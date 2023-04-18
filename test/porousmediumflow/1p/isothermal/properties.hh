// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
 * \ingroup OnePTests
 * \brief A test problem for the one-phase model:
 * water is flowing from bottom to top through and around a low permeable lens.
 */

#ifndef DUMUX_1PTEST_PROBLEM_PROPERTIES_HH
#define DUMUX_1PTEST_PROBLEM_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/1p/model.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "problem.hh"
#include "spatialparams.hh"

#if FORCHHEIMER
#include <dumux/flux/forchheimerslaw.hh>
#endif

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct OnePTest { using InheritsFrom = std::tuple<OneP>; };
struct OnePTestBox { using InheritsFrom = std::tuple<BoxModel, OnePTest>; };
struct OnePTestCCTpfa { using InheritsFrom = std::tuple<CCTpfaModel, OnePTest>; };
struct OnePTestCCMpfa { using InheritsFrom = std::tuple<CCMpfaModel, OnePTest>; };
} // end namespace TTag

// Specialize the fluid system type for this type tag
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePTest>
{
    using Scalar = GetPropType<TypeTag, Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// Specialize the grid type for this type tag
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePTest>
{ using type = Dune::YaspGrid<2>; };

// Specialize the problem type for this type tag
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePTest>
{ using type = OnePTestProblem<TypeTag>; };

// Specialize the spatial params type for this type tag
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePTest>
{
    using GridGeometry = GetPropType<TypeTag, GridGeometry>;
    using Scalar = GetPropType<TypeTag, Scalar>;
    using type = OnePTestSpatialParams<GridGeometry, Scalar>;
};

#ifdef FORCHHEIMER
// Specialize the advection type for this type tag
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::OnePTest>
{ using type = ForchheimersLaw<TypeTag>; };
#endif

} // end namespace Dumux::Properties

#endif
