// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
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
 *
 * \brief test problem for the sequential one-phase model.
 */
#ifndef DUMUX_TEST_1P_PROBLEM_PROPERTIES_HH
#define DUMUX_TEST_1P_PROBLEM_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/porousmediumflow/1p/sequential/diffusion/cellcentered/pressureproperties.hh>
#include <dumux/porousmediumflow/1p/sequential/diffusion/problem.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/velocity.hh>

#include "test_1pproblem.hh"
#include "test_1pspatialparams.hh"

namespace Dumux::Properties
{
// Create new type tags
namespace TTag {
struct TestOneP { using InheritsFrom = std::tuple<FVPressureOneP>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TestOneP> { using type = Dune::YaspGrid<2>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TestOneP>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TestOneP> { using type = TestOnePSpatialParams<TypeTag>; };

//Set the problem
template<class TypeTag>
struct Problem<TypeTag, TTag::TestOneP> { using type = TestProblemOneP<TypeTag>; };
} //end namespace Properties

#endif
