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

class ExternalQuadrupole : public ExternalPotential
    {
    public:
    //! Constructor
    ExternalQuadrupole(std::shared_ptr<SystemDefinition> sysdef,
                        float dg)
        : ExternalPotential(sysdef), m_params(sysdef->getParticleData()->getNTypes())
        {   
            setElectrodeGap(dg);
        }

    //! Destructor
    virtual ~ExternalQuadrupole() { }

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
        
        if (kt != 0.0){
            Scalar d2 = dot(r_i, r_i) / m_electrode_gap / m_electrode_gap;
            energy += Scalar(0.5) * kt * d2;
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
            m_k_translational = v["k"].cast<LongReal>();
        };

        /// Convert a parameter set to a dictionary.
        pybind11::dict asDict()
        {
            pybind11::dict pydict;
            // TODO; pack per-type quantities from the ParamType struct to the Python dictionary.
            pydict["k"] = m_k_translational;
            return pydict;
        };

        LongReal m_k_translational; // strength of the harmonic potentials in kT units
        };

    std::vector<ParamType> m_params;                 // Stored indexed properties (k_translational, k_rotational, m_symmetry) for each particle type
    Scalar m_electrode_gap;                           // Stays global (shared geometry)
    };

namespace detail
    {
void export_ExternalQuadrupole(pybind11::module& m)
    {
    pybind11::class_<ExternalQuadrupole,
                     ExternalPotential,
                     std::shared_ptr<ExternalQuadrupole>>(m, "ExternalQuadrupole")
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
                            float>())
        .def("setParams", &ExternalQuadrupole::setParamsPython)
        .def("getParams", &ExternalQuadrupole::getParamsPython)
        .def_property("electrode_gap",
                      &ExternalQuadrupole::getElectrodeGap,
                      &ExternalQuadrupole::setElectrodeGap);
    }
    } // end namespace detail
    } // namespace hpmc
    } // end namespace hoomd
