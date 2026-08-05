// Copyright (c) 2009-2026 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// TODO: Include the header files of classes that will be exported to Python.

#include <pybind11/pybind11.h>
#include "ExternalAnypole.h"
#include "ExternalCoplane.h"

using namespace hoomd::hpmc::detail;

// TODO: Set the name of the python module to match ${COMPONENT_NAME} (set in
// CMakeLists.txt), prefixed with an underscore.
PYBIND11_MODULE(_dep, m)
    {
        // TODO: Call export_Class(m) for each C++ class to be exported to Python.
        export_ExternalAnypole(m);
        export_ExternalCoplane(m);

#ifdef ENABLE_HIP
        // TODO: Call export_ClassGPU(m) for each GPU enabled C++ class to be exported
        // to Python.
#endif
    }
