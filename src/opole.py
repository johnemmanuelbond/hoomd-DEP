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

@hoomd.logging.modify_namespace(("hpmc", "external"))
class ExternalOctupole(External):
    """
    TODO: Document your component.
    """

    __doc__ = inspect.cleandoc(__doc__).replace(
        "{inherited}", inspect.cleandoc(External._doc_inherited)
    )

    _cpp_class_name = 'ExternalOctupole'
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
                k_perp=object,
                k_para=object,
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
        self._cpp_obj = _dep.ExternalOctupole(
            self._simulation.state._cpp_sys_def,
            self.electrode_gap,
            self.electrode_orientation
        )
        return self._cpp_obj

    def _attach_hook(self):
        # Pass only global setups to the C++ constructor on attachment
        self._cpp_obj = _dep.ExternalOctupole(
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
            raw_k1 = data['k_perp']
            raw_k2 = data['k_para']
            
            k1 = raw_k1(timestep) if isinstance(raw_k1, variant.Variant) else float(raw_k1)
            k2 = raw_k2(timestep) if isinstance(raw_k2, variant.Variant) else float(raw_k2)
            
            # Send the clean numeric values down to C++ vector space
            self._cpp_obj.setParamsCpp(type_id, k1, k2)



class ForceOctupole(Custom):
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
                k_para=object,
                k_perp=object,
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
                raw_k1 = type_data['k_para']
                raw_k2 = type_data['k_perp']
                
                k1 = raw_k1(timestep) if isinstance(raw_k1, variant.Variant) else float(raw_k1)
                k2 = raw_k2(timestep) if isinstance(raw_k2, variant.Variant) else float(raw_k2)

                if k1 !=0 or k2!=0:
                    arrays.force[idx] = ForceOctupole.compute_forces(
                        snap.particles.position[idx], k1, k2, q0, dg
                    )

    @classmethod
    def compute_forces(cls, positions, k_para, k_perp, q0, dg):
        """Compute the forces for a given set of particle positions."""
        cosq, sinq = np.cos(q0), np.sin(q0)
        positions = np.atleast_2d(positions)
        
        # Calculate the projection along the electrode orientation
        f1s = -1.0 * k_para / dg * ( positions[:, 0] * cosq + positions[:, 1] * sinq) / dg
        f2s = -1.0 * k_perp / dg * (-positions[:, 0] * sinq + positions[:, 1] * cosq) / dg

        # Assemble N x 3 output layout
        forces = np.zeros_like(positions)
        forces[:, 0] = (f1s * cosq) + (f2s * -sinq)
        forces[:, 1] = (f1s * sinq) + (f2s *  cosq)
        return forces