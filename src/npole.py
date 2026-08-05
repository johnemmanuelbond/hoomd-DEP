# Copyright (c) 2009-2026 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Harmonic potential that restrains particles to a lattice."""

import hoomd
import numpy as np
import inspect

from hoomd.hpmc.external import External
from hoomd.md.force import Custom

from hoomd.data import TypeParameter
from hoomd.data.parameterdicts import TypeParameterDict
from hoomd import variant

from . import _dep

def _quat_to_angle(quats):
    """Convert an array of quaternions to an array of angles."""
    return 2.0 * np.arctan2(quats[:,-1], quats[:,0])

@hoomd.logging.modify_namespace(("hpmc", "external"))
class ExternalAnypole(External):
    """
    TODO: Document your component.
    """

    __doc__ = inspect.cleandoc(__doc__).replace(
        "{inherited}", inspect.cleandoc(External._doc_inherited)
    )

    _cpp_class_name = 'ExternalAnypole'
    _ext_module = _dep

    def __init__(
        self,
        electrode_gap:float,
        electrode_orientation:variant.variant_like,
    ):
        super().__init__()

        # type dependent configurations
        param_spec = TypeParameter(
            name='params',
            type_kind='particle_types',
            param_dict=TypeParameterDict(
                k_trans=object,
                k_rot=object,
                m_sym=int,
                len_keys=1,
                lenient=True
            )
        )
        self._add_typeparam(param_spec)
        
        # Non-type dependent configurations
        self._q0 = electrode_orientation if isinstance(electrode_orientation, variant.Variant) else variant.Constant(electrode_orientation)
        self._dg = electrode_gap
    
    @property
    def electrode_orientation(self) -> variant.Variant:
        """The global electrode orientation angle (radians)."""
        return self._q0 if isinstance(self._q0, variant.Variant) else variant.Constant(self._q0)

    @electrode_orientation.setter
    def electrode_orientation(self, value: variant.variant_like):
        self._q0 = value

    @property
    def electrode_gap(self) -> float:
        """The global electrode gap distance (length units)."""
        return float(self._dg)

    @electrode_gap.setter
    def electrode_gap(self, value:float):
        self._dg = value

    def _make_cpp_obj(self):
        self._cpp_obj = _dep.ExternalAnypole(
            self._simulation.state._cpp_sys_def,
            self.electrode_gap,
            self.electrode_orientation
        )
        return self._cpp_obj

    def _attach_hook(self):
        # Pass only global setups to the C++ constructor on attachment
        self._cpp_obj = _dep.ExternalAnypole(
            self._simulation.state._cpp_sys_def,
            self.electrode_gap,
            self.electrode_orientation
        )
        super()._attach_hook()

    def _apply_param_update(self, timestep):
        """Resolves variant objects into standard floats right before C++ loop calls."""
        for type_name in self.params:
            type_id = self._simulation.state.particle_types.index(type_name)
            data = self.params[type_name]
            
            # Unpack and dynamically resolve values if they are variant timelines
            raw_kt = data['k_trans']
            raw_kr = data['k_rot']
            
            kt = raw_kt(timestep) if isinstance(raw_kt, variant.Variant) else float(raw_kt)
            kr = raw_kr(timestep) if isinstance(raw_kr, variant.Variant) else float(raw_kr)
            m = int(data['m_sym'])
            
            # Send the clean numeric values down to C++ vector space
            self._cpp_obj.setParamsCpp(type_id, kt, kr, m)



class ForceAnypole(Custom):
    """
    TODO: Document your component.
    """
    def __init__(
        self,
        electrode_gap:float,
        electrode_orientation:variant.variant_like,
    ):
        super().__init__(aniso=True)

        # type dependent configurations
        param_spec = TypeParameter(
            name='params',
            type_kind='particle_types',
            param_dict=TypeParameterDict(
                k_trans=object,
                k_rot=object,
                m_sym=int,
                len_keys=1,
                lenient=True
            )
        )
        self._params = param_spec
        
        # Non-type dependent configurations
        self._q0 = electrode_orientation if isinstance(electrode_orientation, variant.Variant) else variant.Constant(electrode_orientation)
        self._dg = electrode_gap
    
    @property
    def params(self) -> TypeParameter:
        """The type-dependent parameters for the force."""
        return self._params

    @params.setter
    def params(self, value: TypeParameter):
        self._params = value

    @property
    def electrode_orientation(self) -> variant.Variant:
        """The global electrode orientation angle (radians)."""
        return self._q0 if isinstance(self._q0, variant.Variant) else variant.Constant(self._q0)

    @electrode_orientation.setter
    def electrode_orientation(self, value: variant.variant_like):
        self._q0 = value

    @property
    def electrode_gap(self) -> float:
        """The global electrode gap distance (length units)."""
        return float(self._dg)

    @electrode_gap.setter
    def electrode_gap(self, value:float):
        self._dg = value

    def set_forces(self, timestep):
        with self._state.cpu_local_snapshot as snap, self.cpu_local_force_arrays as arrays:
            N = len(snap.particles.position)
            if N == 0:
                return # Safe exit if this MPI rank has no local particles

            # Zero out the output buffers up to N
            arrays.force[:N] = 0.0
            arrays.torque[:N] = 0.0
            arrays.potential_energy[:N] = 0.0

            # Resolve the global electrode angle variant for this exact timestep
            q0 = self.electrode_orientation(timestep)
            dg = self.electrode_gap

            # Loop over all unique system type names to cleanly apply vectorization masks
            for i,t in enumerate(self._state.particle_types):
                idx = snap.particles.typeid == i
                if not any(idx): continue # safe exit if this type has no local particles

                type_data = self.params[t]

                # Safe-unpack variant vs float parameters
                raw_kt = type_data['k_trans']
                raw_kr = type_data['k_rot']
                
                kt = raw_kt(timestep) if isinstance(raw_kt, variant.Variant) else float(raw_kt)
                kr = raw_kr(timestep) if isinstance(raw_kr, variant.Variant) else float(raw_kr)
                m = int(type_data['m_sym'])

                if kt != 0.0:
                    arrays.force[idx]  = ForceAnypole.compute_forces(
                        snap.particles.position[idx], kt, q0, dg
                    )
                if kr != 0.0 and m != 0:
                    arrays.torque[idx] = ForceAnypole.compute_torques(
                        snap.particles.orientation[idx], kr, q0, m
                    )

    @classmethod
    def compute_forces(cls, positions, kt, q0, dg):
        """Compute the forces for a given set of particle positions."""
        cosq, sinq = np.cos(q0), np.sin(q0)
        positions = np.atleast_2d(positions)
        
        # Calculate the projection along the electrode orientation
        drs = (positions[:, 0] * cosq + positions[:, 1] * sinq) / dg
        fs = -1.0 * kt / dg * drs
        
        # Assemble N x 3 output layout
        forces = np.zeros_like(positions)
        forces[:, 0] = fs * cosq
        forces[:, 1] = fs * sinq
        return forces

    @classmethod
    def compute_torques(cls, orientations, kr, q0, m):
        """Compute the torques for a given set of particle orientations (angles)."""
        # Translate internal quaternions  to planar angles relative to electrode orientation
        dts = _quat_to_angle(orientations) - q0

        # Assemble N x 3 torque array layout (torques align strictly along z-axis)
        torques = np.zeros_like(orientations)[:, :3]
        torques[:, 2] = -0.5 * kr * np.sin(2.0 * m * dts)
        return torques

    # @classmethod
    # def add_to_logger(cls, logger, kt, kr, dg):
    #     logger[('log','k_trans')] = (lambda: np.array([kt]), 'scalar')
    #     logger[('log','k_rot')] = (lambda: np.array([kr]), 'scalar')
    #     logger[('log','direct')] = (lambda: np.array([0.0]), 'scalar')
    #     logger[('log','dg')] = (lambda: [dg], 'scalar')