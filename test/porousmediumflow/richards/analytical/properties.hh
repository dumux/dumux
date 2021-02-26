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
 * \ingroup RichardsTests
 * \brief The properties of the one-dimensional infiltration problem with a smooth, given solution.
 */
#ifndef DUMUX_RICHARDS_ANALYTICALPROPERTIES_HH
#define DUMUX_RICHARDS_ANALYTICALPROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "spatialparams.hh"
#include "problem.hh"

//////////
// Specify the properties for the analytical problem
//////////
namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RichardsAnalytical { using InheritsFrom = std::tuple<Richards>; };
struct RichardsAnalyticalBox { using InheritsFrom = std::tuple<RichardsAnalytical, BoxModel>; };
struct RichardsAnalyticalCC { using InheritsFrom = std::tuple<RichardsAnalytical, CCTpfaModel>; };
} // end namespace TTag

// Use 2d YaspGrid
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsAnalytical> { using type = Dune::YaspGrid<2>; };

// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsAnalytical> { using type = RichardsAnalyticalProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsAnalytical>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RichardsAnalyticalSpatialParams<GridGeometry, Scalar>;
};

} // end namespace Dumux::Properties

#endif
