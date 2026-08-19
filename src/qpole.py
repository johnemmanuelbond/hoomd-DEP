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
class ExternalQuadrupole(External):
    """
    TODO: Document your component.
    """

    __doc__ = inspect.cleandoc(__doc__).replace(
        "{inherited}", inspect.cleandoc(External._doc_inherited)
    )

    _cpp_class_name = 'ExternalQuadrupole'
    _ext_module = _dep

    def __init__(
        self,
        electrode_gap:float,
    ):
        super().__init__()

        # type dependent configurations
        param_spec = TypeParameter(
            name='params',
            type_kind='particle_types',
            param_dict=TypeParameterDict(
                k=object,
                len_keys=1,
                lenient=True
            )
        )
        self._add_typeparam(param_spec)
        
        # Non-type dependent configurations
        self._dg = electrode_gap

    @property
    def electrode_gap(self) -> float:
        """The global electrode gap distance (length units)."""
        return float(self._dg)

    @electrode_gap.setter
    def electrode_gap(self, value:float):
        self._dg = value

    def _make_cpp_obj(self):
        self._cpp_obj = _dep.ExternalQuadrupole(
            self._simulation.state._cpp_sys_def,
            self.electrode_gap
        )
        return self._cpp_obj

    def _attach_hook(self):
        # Pass only global setups to the C++ constructor on attachment
        self._cpp_obj = _dep.ExternalQuadrupole(
            self._simulation.state._cpp_sys_def,
            self.electrode_gap
        )
        super()._attach_hook()

    def _apply_param_update(self, timestep):
        """Resolves variant objects into standard floats right before C++ loop calls."""
        for type_name in self.params:
            type_id = self._simulation.state.particle_types.index(type_name)
            data = self.params[type_name]
            
            # Unpack and dynamically resolve values if they are variant timelines
            raw_kt = data['k']
            
            kt = raw_kt(timestep) if isinstance(raw_kt, variant.Variant) else float(raw_kt)
            
            # Send the clean numeric values down to C++ vector space
            self._cpp_obj.setParamsCpp(type_id, kt)



class ForceQuadrupole(Custom):
    """
    TODO: Document your component.
    """
    def __init__(
        self,
        electrode_gap:float,
    ):
        super().__init__(aniso=True)

        # type dependent configurations
        param_spec = TypeParameter(
            name='params',
            type_kind='particle_types',
            param_dict=TypeParameterDict(
                k=object,
                len_keys=1,
                lenient=True
            )
        )
        self._params = param_spec
        
        # Non-type dependent configurations
        self._dg = electrode_gap
    
    @property
    def params(self) -> TypeParameter:
        """The type-dependent parameters for the force."""
        return self._params

    @params.setter
    def params(self, value: TypeParameter):
        self._params = value

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
            dg = self.electrode_gap

            # Loop over all unique system type names to cleanly apply vectorization masks
            for i,t in enumerate(self._state.particle_types):
                idx = snap.particles.typeid == i
                if not any(idx): continue # safe exit if this type has no local particles

                type_data = self.params[t]

                # Safe-unpack variant vs float parameters
                raw_kt = type_data['k']
                kt = raw_kt(timestep) if isinstance(raw_kt, variant.Variant) else float(raw_kt)
                if kt != 0.0:
                    arrays.force[idx]  = ForceQuadrupole.compute_forces(
                        snap.particles.position[idx], kt, dg
                    )

    @classmethod
    def compute_forces(cls, positions, kt, dg):
        """Compute the forces for a given set of particle positions."""
        positions = np.atleast_2d(positions)
        
        # Assemble N x 3 output layout
        forces = np.zeros_like(positions)
        forces[:, 0] = -1.0 * kt / dg * positions[:, 0] / dg
        forces[:, 1] = -1.0 * kt / dg * positions[:, 1] / dg
        return forces