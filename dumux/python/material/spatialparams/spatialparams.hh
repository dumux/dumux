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

#ifndef DUMUX_PYTHON_MATERIAL_SPATIAL_PARAMS_HH
#define DUMUX_PYTHON_MATERIAL_SPATIAL_PARAMS_HH


#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

namespace Dumux::Python {

template <class SP, class... options>
void registerSpatialParams(pybind11::handle scope,
                           pybind11::class_<SP, options...> cls)
{
    using pybind11::operator""_a;
    using GridGeometry = std::decay_t<decltype(std::declval<SP>().gridGeometry())>;

    cls.def(pybind11::init([](std::shared_ptr<GridGeometry> gridGeometry){
        return std::make_shared<SP>(gridGeometry);
    }), "gridGeometry"_a);
}

} // namespace Dumux::Python

#endif
