// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
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
