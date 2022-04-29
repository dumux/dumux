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

#ifndef DUMUX_PYTHON_POROUSMEDIUMFLOW_VELOCITYOUTPUT_HH
#define DUMUX_PYTHON_POROUSMEDIUMFLOW_VELOCITYOUTPUT_HH

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/common/typeregistry.hh>
#include <dumux/porousmediumflow/velocityoutput.hh>

namespace Dumux::Python {

// Python wrapper for the PorousMediumFlowVelocityOutput class
template<class GridVariables, class FluxVariables, class... options>
void registerPorousMediumFlowVelocityOutput(pybind11::handle scope,
                                            pybind11::class_<Dumux::PorousMediumFlowVelocityOutput<GridVariables, FluxVariables>, options...> cls)
{
    using pybind11::operator""_a;
    using namespace Dune::Python;
    using VelocityOutput = Dumux::PorousMediumFlowVelocityOutput<GridVariables, FluxVariables>;

    cls.def(pybind11::init([](const GridVariables& gridVariables){
        return std::make_shared<VelocityOutput>(gridVariables);
    }));
}

} // end namespace Dumux::Python

#endif
