#include <pybind11/pybind11.h>
#include "ExternalQuadrupole.h"
#include "ExternalOctupole.h"
#include "ExternalAnypole.h"
#include "ExternalCoplane.h"

using namespace hoomd::hpmc::detail;

PYBIND11_MODULE(_dep, m)
    {
        export_ExternalQuadrupole(m);
        export_ExternalOctupole(m);
        export_ExternalAnypole(m);
        export_ExternalCoplane(m);

#ifdef ENABLE_HIP
        // TODO: Call export_ClassGPU(m) for each GPU enabled C++ class to be exported
        // to Python.
#endif
    }
