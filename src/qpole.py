# -*- coding: utf-8 -*-
"""
Contains external potentials (for HPMC) and external forces (for MD) that model
dielectrophoresis (DEP) of colloidal particles in a quadrupolar electrode
geometry. The effective interaction can be represented as a harmonic trap within
the plane of the electrodes:

.. math::

    U = \\frac{1}{2} k r^2 / d_g^2

Where :math:`d_g` is the electrode gap and :math:`k` parameterizes the energy scale strength
of the particle field interaction, in :math:`k_BT` inits. Most generally, :math:`k` depends
on the particle geometry, electrode geometry, the material properties of the particle and medium, 
and the applied field strength. Most practically, it is proportional to :math:`V_{pp}^2`, where 
:math:`V_{pp}` is the AC voltage drop across the electrode. Since :math:`k` depends on particle 
properties as well as electrode properties, both :py:class:`ExternalQuadrupole` and 
:py:class:`ForceQuadrupole` allow :math:`k` to be specified on a per-particle-type basis. 
"""

import hoomd
import numpy as np
import inspect

from hoomd.hpmc.external import External
from hoomd.md.force import Custom

from hoomd.data import TypeParameter
from hoomd.data.parameterdicts import TypeParameterDict
from hoomd import variant

try: from . import _dep
except ImportError:
    from unittest.mock import MagicMock
    _dep = MagicMock()

@hoomd.logging.modify_namespace(("hpmc", "external", "ExternalQuadrupole"))
class ExternalQuadrupole(External):
    """
    Apply a harmonic potential energy in the xy plane, centered at the origin.

    Args:
        electrode_gap (float): Separation between the electrodes in length units.

    The field is defined by a single global electrode spacing, ``electrode_gap``,
    and a per-type interaction strength ``params[type].k``. The interaction acts
    only in the in-plane directions, with zero contribution along the z axis.

    Example:
        The type parameter is specified as a dictionary keyed by particle type:

        .. code-block:: python

            external = hoomd.dep.ExternalQuadrupole(electrode_gap=100.0)
            external.params["A"] = {"k": 250}
    
    Note:
        `ExternalQuadrupole` does not support execution on GPUs.

    {inherited}
        
    **Members defined in** `ExternalQuadrupole`

    Attributes:
        electrode_gap (float): Separation between the electrodes in length units.

    .. py:attribute:: params

        Per-particle-type quadrupole coefficients. The dictionary has the following keys:

        * ``k``: (`float` or `variant-like`, **required**) - :math:`k` :math:`[\\mathrm{energy}]`

        Type: `TypeParameter` [``particle_type``, `dict`]
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
    Apply a harmonic restoring force in the xy plane, centered at the origin.

    Args:
        electrode_gap (float): Separation between the electrodes in length units.

    The field is defined by a single global electrode spacing, ``electrode_gap``,
    and a per-type interaction strength ``params[type].k``. The force is applied
    only in the in-plane directions, with zero force along the z axis.

    Example:
        .. code-block:: python

            force = hoomd.dep.ForceQuadrupole(electrode_gap=100.0)
            force.params["A"] = {"k": 250}
    
    {inherited}

    **Members defined in** `ForceQuadrupole`
            
    Attributes:
        electrode_gap (float): Separation between the electrodes in length units.

    .. py:attribute:: params

        Per-particle-type quadrupole coefficients. The dictionary has the following keys:

        * ``k``: (`float` or `variant-like`, **required**) - :math:`k` :math:`[\\mathrm{energy}]`

        Type: `TypeParameter` [``particle_type``, `dict`]

    """

    __doc__ = inspect.cleandoc(__doc__).replace(
        "{inherited}", inspect.cleandoc(Custom._doc_inherited)
    )

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
        return self._params

    @params.setter
    def params(self, value: TypeParameter):
        self._params = value

    @property
    def electrode_gap(self) -> float:
        return float(self._dg)

    @electrode_gap.setter
    def electrode_gap(self, value:float):
        self._dg = value

    def set_forces(self, timestep):
        """Set the forces in the simulation loop.

        Args:
            timestep (int): The current timestep in the simulation.
        """
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
        """
        Compute the dielectrophoretic forces on each particle given their positions.

        Args:
            positions (np.ndarray): An N x 3 array of particle positions.
            kt (float): The per-type quadrupole strength.
            dg (float): The global electrode gap distance.
        Returns:
            np.ndarray: An N x 3 array of forces.
        """
        positions = np.atleast_2d(positions)
        
        # Assemble N x 3 output layout
        forces = np.zeros_like(positions)
        forces[:, 0] = -1.0 * kt / dg * positions[:, 0] / dg
        forces[:, 1] = -1.0 * kt / dg * positions[:, 1] / dg
        return forces