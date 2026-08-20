# -*- coding: utf-8 -*-
"""
Contains external potentials (for HPMC) and external forces (for MD) that model
dielectrophoresis (DEP) of colloidal particles in an octupolar electrode
geometry. The effective interaction can be approximated as a pair of orthogonal harmonic
traps within the plane of the electrodes, with an additional orientational dependence
set by the electrode angle:

.. math::

    U/k_BT = \\frac{1}{2} k_{\\parallel} (\\mathbf{r}\\cdot\\hat{\\mathbf{e}})^2/d_g^2 +
    \\frac{1}{2} k_{\\perp} (\\mathbf{r}\\times\\hat{\\mathbf{e}})^2 / d_g^2

Where :math:`d_g` is the electrode gap, :math:`\\hat{\\mathbf{e}}` is the
orientation of the electrode field, and the coefficients :math:`k_{\\parallel}`
and :math:`k_{\\perp}` parameterize the energy scale, in :math:`k_BT` units, of 
the particle-field interaction along and transverse to the electrode orientation. 
Most generally, these coefficients depend on the particle geometry, electrode geometry, 
the material properties of the particle and medium, and the applied field strength.
Practically, they are proportional to :math:`V_{i}^2`, where :math:`V_{\\parallel}`
and :math:`V_{\\perp}` are the AC voltage drops across the two electrode channels.
Since these coefficients depend on particle properties as well as electrode properties,
both :py:class:`ExternalOctupole` and :py:class:`ForceOctupole` allow
:math:`k_{\\parallel}` and :math:`k_{\\perp}` to be specified on a
per-particle-type basis.
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

@hoomd.logging.modify_namespace(("hpmc", "external", "ExternalOctupole"))
class ExternalOctupole(External):
    """
    Apply an orientable octupolar field between a pair of electrodes.

    Args:
        electrode_gap (float): Separation between the electrodes in length units.
        electrode_orientation (float or hoomd.variant.Variant): Global electrode
            orientation angle in radians.

    The interaction is constructed from a global gap, ``electrode_gap``, a
    global orientation angle, ``electrode_orientation``, and per-type interaction
    strengths ``params[type].k_para`` and ``params[type].k_perp``.

    Example:
        .. code-block:: python

            external = hoomd.dep.ExternalOctupole(
                electrode_gap=100.0,
                electrode_orientation=3.1415/4,
            )
            external.params["A"] = {"k_para": 250, "k_perp": 100}

    Note:
        `ExternalOctupole` does not support execution on GPUs.

    {inherited}

    **Members defined in** `ExternalOctupole`

    Attributes:
        electrode_gap (float): Separation between the electrodes in length units.
        electrode_orientation (hoomd.variant.Variant): Global electrode angle.

    .. py:attribute:: params
    
        Per-particle-type octupole coefficients. The dictionary has the following keys:

        * ``k_para``: (`float` or `variant-like`, **required**) - :math:`k_{\\parallel}` :math:`[\\mathrm{energy}]`
        * ``k_perp``: (`float` or `variant-like`, **required**) - :math:`k_{\\perp}` :math:`[\\mathrm{energy}]`

        Type: `TypeParameter` [``particle_type``, `dict`]
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
        return self._q0 if isinstance(self._q0, variant.Variant) else variant.Constant(self._q0)

    @electrode_orientation.setter
    def electrode_orientation(self, value: variant.variant_like):
        self._q0 = value

    @property
    def electrode_gap(self) -> float:
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
    Apply an orientable octupolar force in the plane of the electrodes.

    Args:
        electrode_gap (float): Separation between the electrodes in length units.
        electrode_orientation (float or hoomd.variant.Variant): Global electrode
            orientation angle in radians.

    The interaction is constructed from a global gap, ``electrode_gap``, a
    global orientation angle, ``electrode_orientation``, and per-type interaction
    strengths ``params[type].k_para`` and ``params[type].k_perp``.

    Example:
        .. code-block:: python

            force = hoomd.dep.ForceOctupole(
                electrode_gap=100.0,
                electrode_orientation=3.1415/4,
            )
            force.params["A"] = {"k_para": 250, "k_perp": 100}

    {inherited}

    Attributes:
        electrode_gap (float): Separation between the electrodes in length units.
        electrode_orientation (hoomd.variant.Variant): Global electrode angle.
    
    .. py:attribute:: params
    
        Per-particle-type octupole coefficients. The dictionary has the following keys:

        * ``k_para``: (`float` or `variant-like`, **required**) - :math:`k_{\\parallel}` :math:`[\\mathrm{energy}]`
        * ``k_perp``: (`float` or `variant-like`, **required**) - :math:`k_{\\perp}` :math:`[\\mathrm{energy}]`

        Type: `TypeParameter` [``particle_type``, `dict`]
    """

    __doc__ = inspect.cleandoc(__doc__).replace(
        "{inherited}", inspect.cleandoc(Custom._doc_inherited)
    )

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
        return self._params

    @params.setter
    def params(self, value: TypeParameter):
        self._params = value

    @property
    def electrode_orientation(self) -> variant.Variant:
        return self._q0 if isinstance(self._q0, variant.Variant) else variant.Constant(self._q0)

    @electrode_orientation.setter
    def electrode_orientation(self, value: variant.variant_like):
        self._q0 = value

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
        """Compute the dielectrophoretic forces on each particle given their positions.

        Args:
            positions (np.ndarray): An N x 3 array of particle positions.
            k_para (float): Strength of the force component parallel to the
                electrode orientation.
            k_perp (float): Strength of the force component perpendicular to the
                electrode orientation.
            q0 (float): Global electrode orientation angle in radians.
            dg (float): The global electrode gap distance.

        Returns:
            np.ndarray: An N x 3 array of forces.
        """
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