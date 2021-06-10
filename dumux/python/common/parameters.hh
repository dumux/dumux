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

#ifndef DUMUX_PYTHON_COMMON_PARAMETERS_HH
#define DUMUX_PYTHON_COMMON_PARAMETERS_HH

#include <dumux/common/parameters.hh>

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

namespace Dumux::Python {

template <class... options>
void registerParameters(pybind11::handle scope,
                        pybind11::class_<Parameters, options...> cls)
{
    cls.def(pybind11::init());

    cls.def("init", [](Parameters& self) { self.init(); });
    cls.def("init", [](Parameters& self, const std::string& parameterFileName) { self.init(parameterFileName); });
}

template<class Scalar>
void registerParameters(pybind11::handle scope, const char *clsName = "Parameters")
{
    pybind11::class_<Parameters> cls(scope, clsName);
    registerParameters(scope, cls);
}

} // namespace Dumux::Python

#endif
