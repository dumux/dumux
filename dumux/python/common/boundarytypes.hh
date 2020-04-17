#ifndef DUMUX_PYTHON_COMMON_BOUNDARYTYPES_HH
#define DUMUX_PYTHON_COMMON_BOUNDARYTYPES_HH

#include <dune/common/classname.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/common/typeregistry.hh>

#include <dumux/common/boundarytypes.hh>

namespace Dumux::Python {

template <class BoundaryTypes, class... Options>
void registerBoundaryTypes(pybind11::handle scope, pybind11::class_<BoundaryTypes, Options...> cls)
{
    using pybind11::operator""_a;

    cls.def(pybind11::init());
    cls.def("reset", &BoundaryTypes::reset);
    cls.def("setNeumann", &BoundaryTypes::setAllNeumann);
    cls.def("setDirichlet", &BoundaryTypes::setAllDirichlet);
    cls.def("isDirichlet", &BoundaryTypes::hasDirichlet);
    cls.def("isNeumann", &BoundaryTypes::hasNeumann);
}

template <class BoundaryTypes>
void registerBoundaryTypes(pybind11::handle scope)
{
    using namespace Dune::Python;

    auto [cls, addedToRegistry] = insertClass<BoundaryTypes>(
        scope, "BoundaryTypes",
        GenerateTypeName(Dune::className<BoundaryTypes>()),
        IncludeFiles{"dumux/python/common/boundarytypes.hh"}
    );

    if (addedToRegistry)
        registerBoundaryTypes(scope, cls);
}

} // namespace Dumux::Python

#endif
