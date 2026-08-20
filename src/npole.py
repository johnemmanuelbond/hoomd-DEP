# -*- coding: utf-8 -*-
"""
Contains external potentials (for HPMC) and external forces (for MD) that model
an abstract multipolar DEP interaction in a generalized electrode geometry. This
module is not intended to reproduce a specific colloidal experiment or a
single known device architecture. Instead, it is designed as a flexible,
future-facing framework that can represent a broad class of unknown or
emergent multipolar systems through an orientable, symmetry-controlled
interaction field.

The effective potential energy can be written as a generalized translational
bias along a field direction together with a symmetry-dependent rotational
coupling:

.. math::

    U = \\frac{1}{2} k_{t} (\\mathbf{r}\\cdot\\hat{\\mathbf{e}})^2/d_g^2
    + \\frac{1}{2} k_{r} \\sin^2(2m(\\theta-\\theta_0))

where :math:`d_g` is the characteristic gap scale, :math:`\\hat{\\mathbf{e}}` is
an orienting axis along :math:`\\theta_0`, :math:`\\theta` is the in-plane orientation 
angle, and :math:`m` is the symmetry order. The coefficients :math:`k_{t}` and
:math:`k_{r}` are treated as effective coupling strengths that can be
adapted to different particle geometries or field configurations. As a result,
both :py:class:`ExternalAnypole` and :py:class:`ForceAnypole` allow these
coefficients to be specified on a per-particle-type basis.
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

@hoomd.logging.modify_namespace(("hpmc", "external", "ExternalAnypole"))
class ExternalAnypole(External):
    """
    Apply a generalized multipolar external potential with a translational bias and
    a symmetry-controlled orientational coupling.

    Args:
        electrode_gap (float): Characteristic separation scale in length units.
        electrode_orientation (float or hoomd.variant.Variant): Global orientation
            angle in radians.

    This class is intentionally abstract: it is not written to reproduce one
    specific experimental architecture, but to provide a flexible template for
    future systems whose effective interactions are described by a generalized
    orientable field. The field is defined by a global gap,
    ``electrode_gap``, a global orientation, ``electrode_orientation``, and
    per-type parameters ``params[type].k_trans``, ``params[type].k_rot``, and
    ``params[type].m_sym``.

    Example:
        .. code-block:: python

            external = hoomd.dep.ExternalAnypole(
                electrode_gap=100.0,
                electrode_orientation=3.1415/4,
            )
            external.params["A"] = {"k_trans": 250, "k_rot": 100, "m_sym": 2}

    Note:
        `ExternalAnypole` does not support execution on GPUs.

    {inherited}

    **Members defined in** `ExternalAnypole`

    Attributes:
        electrode_gap (float): Characteristic separation scale in length units.
        electrode_orientation (hoomd.variant.Variant): Global field orientation.

    .. py:attribute:: params

        Per-particle-type generalized anypole coefficients. The dictionary has the following keys:

        * ``k_trans``: (`float` or `variant-like`, **required**) - :math:`k_t` :math:`[\\mathrm{energy}]`
        * ``k_rot``: (`float` or `variant-like`, **required**) - :math:`k_r` :math:`[\\mathrm{energy}]`
        * ``m_sym``: (`int`, **required**) - :math:`m`

        Type: `TypeParameter` [``particle_type``, `dict`]

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
    Apply a generalized multipolar force with an orientational torque term.

    Args:
        electrode_gap (float): Characteristic separation scale in length units.
        electrode_orientation (float or hoomd.variant.Variant): Global orientation
            angle in radians.

    This force is designed as a general-purpose template for unknown or future
    multipolar systems rather than a strict model of a specific colloidal setup.
    It combines a translational field along the local electrode direction with a
    torque that depends on the particle orientation relative to the field.

    Example:
        .. code-block:: python

            force = hoomd.dep.ForceAnypole(
                electrode_gap=100.0,
                electrode_orientation=3.1415/4,
            )
            force.params["A"] = {"k_trans": 250, "k_rot": 100, "m_sym": 2}

    {inherited}

    **Members defined in** `ForceAnypole`

    Attributes:
        electrode_gap (float): Characteristic separation scale in length units.
        electrode_orientation (hoomd.variant.Variant): Global field orientation.

    .. py:attribute:: params

        Per-particle-type generalized anypole coefficients. The dictionary has the following keys:

        * ``k_trans``: (`float` or `variant-like`, **required**) - :math:`k_t` :math:`[\\mathrm{energy}]`
        * ``k_rot``: (`float` or `variant-like`, **required**) - :math:`k_r` :math:`[\\mathrm{energy}]`
        * ``m_sym``: (`int`, **required**) - :math:`m`

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
        """Compute the generalized translational force along the field direction.

        Args:
            positions (np.ndarray): An N x 3 array of particle positions.
            kt (float): Translational coupling strength.
            q0 (float): Global orientation angle in radians.
            dg (float): The characteristic electrode gap scale.

        Returns:
            np.ndarray: An N x 3 array of forces.
        """
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
        """Compute the symmetry-controlled orientational torque.

        Args:
            orientations (np.ndarray): An N x 4 array of particle quaternions.
            kr (float): Rotational coupling strength.
            q0 (float): Global field orientation angle in radians.
            m (int): Angular symmetry order of the generalized field.

        Returns:
            np.ndarray: An N x 3 array of torques, with the z component active.
        """
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