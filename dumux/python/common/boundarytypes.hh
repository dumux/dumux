// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief TODO: docme!
 */

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
    cls.def("__copy__",  [](const BoundaryTypes& self) {
        return BoundaryTypes(self);
    });
    cls.def("__deepcopy__", [](const BoundaryTypes& self, pybind11::dict) {
        return BoundaryTypes(self);
    }, "memo"_a);
    cls.def("reset", &BoundaryTypes::reset);
    cls.def("setNeumann", &BoundaryTypes::setAllNeumann);
    cls.def("setDirichlet", &BoundaryTypes::setAllDirichlet);
    cls.def_property_readonly("isDirichlet", &BoundaryTypes::hasDirichlet);
    cls.def_property_readonly("isNeumann", &BoundaryTypes::hasNeumann);
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
