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
 * \ingroup TwoPTwoCTests
 * \brief The properties of the non-isothermal gas injection problem where a gas (e.g. air) is injected into a fully water saturated medium.
 */

#ifndef DUMUX_WATER_AIR_PROPERTIES_HH
#define DUMUX_WATER_AIR_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/porousmediumflow/2p2c/model.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct WaterAir { using InheritsFrom = std::tuple<TwoPTwoCNI>; };
struct WaterAirBox { using InheritsFrom = std::tuple<WaterAir, BoxModel>; };
struct WaterAirCCTpfa { using InheritsFrom = std::tuple<WaterAir, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::WaterAir> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::WaterAir> { using type = WaterAirProblem<TypeTag>; };

// Set the wetting phase
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::WaterAir> { using type = FluidSystems::H2ON2<GetPropType<TypeTag, Properties::Scalar>>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::WaterAir>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = WaterAirSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::WaterAir> { static constexpr bool value = true; };

// Enable caching
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::WaterAir> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::WaterAir> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::WaterAir> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
