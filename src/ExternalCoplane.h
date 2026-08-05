// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#pragma once

#include "hoomd/Compute.h"
#include "hoomd/HOOMDMath.h"
#include "hoomd/VectorMath.h"
#include <hoomd/Variant.h>

#include "hoomd/hpmc/ExternalPotential.h"

#ifndef __HIPCC__
#include <pybind11/pybind11.h>
#endif

namespace hoomd
    {
namespace hpmc
    {

class ExternalCoplane : public ExternalPotential
    {
    public:
    //! Constructor
    ExternalCoplane(std::shared_ptr<SystemDefinition> sysdef,
                        float dg)
        : ExternalPotential(sysdef), m_params(sysdef->getParticleData()->getNTypes())
        {   
            setElectrodeGap(dg);
        }

    //! Destructor
    virtual ~ExternalCoplane() { }

    //! Setter for electrode gap
    void setElectrodeGap(const float& dg) { m_electrode_gap = Scalar(dg); }
    float getElectrodeGap() { return static_cast<float>(m_electrode_gap); }
    
    //! Set parameters for a specific particle type
    void setParamsPython(const std::string& particle_type,
                                  pybind11::dict params)
        {
        unsigned int particle_type_id = m_sysdef->getParticleData()->getTypeByName(particle_type);
        m_params[particle_type_id] = ParamType(params);
        }

    //! Get parameters for a specific particle type
    pybind11::dict getParamsPython(const std::string& particle_type)
        {
        unsigned int particle_type_id = m_sysdef->getParticleData()->getTypeByName(particle_type);
        return m_params[particle_type_id].asDict();
        }


    protected:
    virtual LongReal particleEnergyImplementation(uint64_t timestep,
                                                  unsigned int tag_i,
                                                  unsigned int type_i,
                                                  const vec3<LongReal>& r_i,
                                                  const quat<LongReal>& q_i,
                                                  LongReal charge_i,
                                                  Trial trial)
        {
        Scalar energy = 0.0;

        auto& params = m_params[type_i];
        Scalar kt = params.m_k_translational;
        Scalar kr = params.m_k_rotational;
        Scalar m = params.m_m_symmetry;
        
        if (kt != 0.0){
            Scalar d = r_i.x / m_electrode_gap;
            Scalar arg = std::clamp(Scalar(2.45) * d, Scalar(-0.99*M_PI/2), Scalar(0.99*M_PI/2));
            energy += Scalar(0.5) * kt * tan(arg) * tan(arg);
        }

        if (kr != 0.0 && m != 0.0){
            vec3<Scalar> xhat = vec3<Scalar>(1.0, 0.0, 0.0);
            Scalar c = q_i.s;
            Scalar s = fast::sqrt(dot(q_i.v, q_i.v));
            vec3<Scalar> major_axis = xhat + 2*c*cross(q_i.v, xhat) + 2*s*cross(q_i.v, cross(q_i.v, xhat));
            Scalar dtheta = atan2(major_axis.y, major_axis.x);
            energy += Scalar(0.5) * kr * sin(m * dtheta) * sin(m * dtheta);
        }

        return energy;
        }
    
    //! Define a flat plain data structure for per-type parameters
    struct ParamType
        {
        ParamType() { }

        /// Construct a parameter set from a dictionary.
        ParamType(pybind11::dict params){
            pybind11::dict v = params;
            m_k_translational = v["k_trans"].cast<LongReal>();
            m_k_rotational    = v["k_rot"].cast<LongReal>();
            m_m_symmetry      = v["m_sym"].cast<int>();
        };

        /// Convert a parameter set to a dictionary.
        pybind11::dict asDict()
        {
            pybind11::dict pydict;
            // TODO; pack per-type quantities from the ParamType struct to the Python dictionary.
            pydict["k_trans"] = m_k_translational;
            pydict["k_rot"] = m_k_rotational;
            pydict["m_sym"] = m_m_symmetry;
            return pydict;
        };

        LongReal m_k_translational; // strength of the harmonic potentials in kT units
        LongReal m_k_rotational;
        int m_m_symmetry; // number of symmetrically equivalent orientations
        };

    std::vector<ParamType> m_params;                 // Stored indexed properties (k_translational, k_rotational, m_symmetry) for each particle type
    Scalar m_electrode_gap;                           // Stays global (shared geometry)
    };

namespace detail
    {
void export_ExternalCoplane(pybind11::module& m)
    {
    pybind11::class_<ExternalCoplane,
                     ExternalPotential,
                     std::shared_ptr<ExternalCoplane>>(m, "ExternalCoplane")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            float>())
        .def("setParams", &ExternalCoplane::setParamsPython)
        .def("getParams", &ExternalCoplane::getParamsPython)
        .def_property("electrode_gap",
                      &ExternalCoplane::getElectrodeGap,
                      &ExternalCoplane::setElectrodeGap);
    }
    } // end namespace detail
    } // namespace hpmc
    } // end namespace hoomd
