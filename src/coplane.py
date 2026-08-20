# -*- coding: utf-8 -*-
"""
Contains external potentials (for HPMC) and external forces (for MD) that model
dielectrophoresis (DEP) of colloidal particles in a semi-infinite coplanar electrode
geometry. In this configuration, the translational potential energy is not purely
harmonic because the electric field itself rapidly increases near the electrode edges.
Following |simple-models|_, the translational potential energy is well-approximated by:

.. math::

    U/k_BT = \\frac{1}{2} k_{t} \\tan^2(\\beta x / d_g) + \\frac{1}{2} k_{r} \\sin^2(2 m \\theta)

where :math:`d_g` is the electrode gap, :math:`x` is the particle position along
the x-axis, and :math:`\\beta=2.45` is a geometric prefactor. The interaction
also includes a rotational contribution that depends on the in-plane particle
orientation :math:`\\theta` and the symmetry order :math:`m`.

As with the other DEP components, the energy scales :math:`k_{t}`
and :math:`k_{r}`, express in :math:`k_BT` units, depends on particle 
and electrode properties, as well as the applied field strength, and can therefore
be assigned on a per-particle-type basis. Both :py:class:`ExternalCoplane` and
:py:class:`ForceCoplane` expose these parameters through the ``params[type]`` dictionary.

.. |simple-models| replace:: (Zhang, et. al. *Langmuir*, 2024)
.. _simple-models: https://doi.org/10.1021/acs.langmuir.4c03101
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

def _quat_to_angle(quats):
    """Convert an array of quaternions to an array of angles."""
    return 2.0 * np.arctan2(quats[:,-1], quats[:,0])

@hoomd.logging.modify_namespace(("hpmc", "external", "ExternalCoplane"))
class ExternalCoplane(External):
    """
    Apply a coplanar DEP potential that combines translational and rotational
    steering in the plane of the electrodes.

    Args:
        electrode_gap (float): Separation between the electrodes in length units.

    The interaction is constructed from a global electrode gap, ``electrode_gap``,
    per-type interaction strengths ``params[type].k_trans`` and ``params[type].k_rot``,
    and per-type particle symmetries, ``params[type].m_sym``, which denote the number
    of symmetrically equivalent orientations of each particle type.

    Example:
        .. code-block:: python

            external = hoomd.dep.ExternalCoplane(electrode_gap=100.0)
            external.params["A"] = {"k_trans": 250, "k_rot": 100, "m_sym": 1}

    Note:
        `ExternalCoplane` does not support execution on GPUs.

    {inherited}

    **Members defined in** `ExternalCoplane`

    Attributes:
        electrode_gap (float): Separation between the electrodes in length units.

    .. py:attribute:: params

        Per-particle-type generalized coplane coefficients. The dictionary has the following keys:

        * ``k_trans``: (`float` or `variant-like`, **required**) - :math:`k_{t}` :math:`[\\mathrm{energy}]` 
        * ``k_rot``: (`float` or `variant-like`, **required**) - :math:`k_{r}` :math:`[\\mathrm{energy}]` 
        * ``m_sym``: (`int`, **required**) - :math:`m`

        Type: `TypeParameter` [``particle_type``, `dict`]
        
    """

    __doc__ = inspect.cleandoc(__doc__).replace(
        "{inherited}", inspect.cleandoc(External._doc_inherited)
    )

    _cpp_class_name = 'ExternalCoplane'
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
                k_trans=object,
                k_rot=object,
                m_sym=int,
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
        self._cpp_obj = _dep.ExternalCoplane(
            self._simulation.state._cpp_sys_def,
            self.electrode_gap,
        )
        return self._cpp_obj

    def _attach_hook(self):
        # Pass only global setups to the C++ constructor on attachment
        self._cpp_obj = _dep.ExternalCoplane(
            self._simulation.state._cpp_sys_def,
            self.electrode_gap,
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



class ForceCoplane(Custom):
    """
    Apply a coplanar DEP force with an orientational torque term.

    Args:
        electrode_gap (float): Separation between the electrodes in length units.

    The interaction is constructed from a global electrode gap, ``electrode_gap``,
    per-type interaction strengths ``params[type].k_trans`` and ``params[type].k_rot``,
    and per-type particle symmetries, ``params[type].m_sym``, which denote the number
    of symmetrically equivalent orientations of each particle type.

    Example:
        .. code-block:: python

            force = hoomd.dep.ForceCoplane(electrode_gap=100.0)
            force.params["A"] = {"k_trans": 250, "k_rot": 100, "m_sym": 1}

    {inherited}

    Attributes:
        electrode_gap (float): Separation between the electrodes in length units.

    .. py:attribute:: params

        Per-particle-type generalized coplane coefficients. The dictionary has the following keys:


        * ``k_trans``: (`float` or `variant-like`, **required**) - :math:`k_{t}` :math:`[\\mathrm{energy}]`
        * ``k_rot``: (`float` or `variant-like`, **required**) - :math:`k_{r}` :math:`[\\mathrm{energy}]`
        * ``m_sym``: (`int`, **required**) - :math:`m`

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
                k_trans=object,
                k_rot=object,
                m_sym=int,
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
        """Set the forces and torques in the simulation loop.

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
                raw_kt = type_data['k_trans']
                raw_kr = type_data['k_rot']
                
                kt = raw_kt(timestep) if isinstance(raw_kt, variant.Variant) else float(raw_kt)
                kr = raw_kr(timestep) if isinstance(raw_kr, variant.Variant) else float(raw_kr)
                m = int(type_data['m_sym'])

                if kt != 0.0:
                    arrays.force[idx]  = ForceCoplane.compute_forces(
                        snap.particles.position[idx], kt, dg
                    )
                if kr != 0.0 and m != 0:
                    arrays.torque[idx] = ForceCoplane.compute_torques(
                        snap.particles.orientation[idx], kr, m
                    )

    @classmethod
    def compute_forces(cls, positions, kt, dg):
        """Compute the in-plane coplanar DEP force on each particle.

        Args:
            positions (np.ndarray): An N x 3 array of particle positions.
            kt (float): Translational coupling strength.
            dg (float): The global electrode gap distance.

        Returns:
            np.ndarray: An N x 3 array of forces.
        """
        xs = np.array(positions)[:,0]

        beta = 2.45
        arg = np.clip(beta * xs / dg, -0.99*np.pi/2, 0.99*np.pi/2) # Avoid singularities
        fx = -1.0 * kt * beta / dg * np.tan(arg) * (1 / np.cos(arg))**2
        
        # Assemble N x 3 output layout
        forces = np.zeros_like(positions)
        forces[:, 0] = fx
        return forces

    @classmethod
    def compute_torques(cls, orientations, kr, m):
        """Compute the out-of-plane torque from the orientational symmetry term.

        Args:
            orientations (np.ndarray): An N x 4 array of particle quaternions.
            kr (float): Rotational coupling strength.
            m (int): Angular symmetry order of the particles represented by ``orientations``.

        Returns:
            np.ndarray: An N x 3 array of torques, with the z component active.
        """
        # Translate internal quaternions  to planar angles relative to electrode orientation
        dts = _quat_to_angle(orientations)
        # Assemble N x 3 torque array layout (torques align strictly along z-axis)
        torques = np.zeros_like(orientations)[:, :3]
        torques[:, 2] = -0.5 * kr * np.sin(2.0 * m * dts)
        return torques