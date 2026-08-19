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

class ExternalOctupole : public ExternalPotential
    {
    public:
    //! Constructor
    ExternalOctupole(std::shared_ptr<SystemDefinition> sysdef,
                        float dg, std::shared_ptr<Variant> q0)
        : ExternalPotential(sysdef), m_params(sysdef->getParticleData()->getNTypes())
        {   
            setElectrodeOrientation(q0);
            setElectrodeGap(dg);
        }

    //! Destructor
    virtual ~ExternalOctupole() { }

    //! Setter for electrode orientation
    void setElectrodeOrientation(const std::shared_ptr<Variant>& q0) { m_electrode_orientation = q0; }
    std::shared_ptr<Variant> getElectrodeOrientation() { return m_electrode_orientation; }

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
        Scalar kt1 = params.m_k_translational_para;
        Scalar kt2 = params.m_k_translational_perp;
        Scalar q0 = (*m_electrode_orientation)(timestep);
        
        if (kt1 != 0.0){
            vec3<Scalar> uvec = vec3<Scalar>(cos(q0), sin(q0), 0.0);
            Scalar d = dot(r_i, uvec) / m_electrode_gap;
            energy += Scalar(0.5) * kt1 * d * d;
        }

        if (kt2 != 0.0){
            vec3<Scalar> uvec = vec3<Scalar>(sin(q0), -cos(q0), 0.0);
            Scalar d = dot(r_i, uvec) / m_electrode_gap;
            energy += Scalar(0.5) * kt2 * d * d;
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
            m_k_translational_para = v["k_para"].cast<LongReal>();
            m_k_translational_perp = v["k_perp"].cast<LongReal>();
        };

        /// Convert a parameter set to a dictionary.
        pybind11::dict asDict()
        {
            pybind11::dict pydict;
            // TODO; pack per-type quantities from the ParamType struct to the Python dictionary.
            pydict["k_para"] = m_k_translational_para;
            pydict["k_perp"] = m_k_translational_perp;
            return pydict;
        };

        LongReal m_k_translational_para; // strength of the harmonic potentials in kT units
        LongReal m_k_translational_perp;
        };

    std::vector<ParamType> m_params;                 // Stored indexed properties (k_translational, k_rotational, m_symmetry) for each particle type
    std::shared_ptr<Variant> m_electrode_orientation; // Stays global (shared field orientation)
    Scalar m_electrode_gap;                           // Stays global (shared geometry)
    };

namespace detail
    {
void export_ExternalOctupole(pybind11::module& m)
    {
    pybind11::class_<ExternalOctupole,
                     ExternalPotential,
                     std::shared_ptr<ExternalOctupole>>(m, "ExternalOctupole")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            float, std::shared_ptr<Variant>>())
        .def("setParams", &ExternalOctupole::setParamsPython)
        .def("getParams", &ExternalOctupole::getParamsPython)
        .def_property("electrode_gap",
                      &ExternalOctupole::getElectrodeGap,
                      &ExternalOctupole::setElectrodeGap)
        .def_property("electrode_orientation",
                      &ExternalOctupole::getElectrodeOrientation,
                      &ExternalOctupole::setElectrodeOrientation);
    }
    } // end namespace detail
    } // namespace hpmc
    } // end namespace hoomd
