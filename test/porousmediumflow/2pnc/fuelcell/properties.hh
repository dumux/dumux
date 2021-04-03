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
 * \ingroup TwoPNCTests
 * \brief The properties of a problem for water management in PEM fuel cells.
 */

#ifndef DUMUX_FUELCELL_PROPERTIES_HH
#define DUMUX_FUELCELL_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/2pnc/model.hh>
#include <dumux/material/fluidsystems/h2on2o2.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
#ifdef NONISOTHERMAL
struct FuelCell { using InheritsFrom = std::tuple<TwoPNCNI>; };
struct FuelCellNIBox { using InheritsFrom = std::tuple<FuelCell, BoxModel>; };
#else
struct FuelCell { using InheritsFrom = std::tuple<TwoPNC>; };
struct FuelCellBox { using InheritsFrom = std::tuple<FuelCell, BoxModel>; };
struct FuelCellCCTpfa { using InheritsFrom = std::tuple<FuelCell, CCTpfaModel>; };
#endif
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::FuelCell> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::FuelCell> { using type = FuelCellProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::FuelCell>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FuelCellSpatialParams<GridGeometry, Scalar>;
};

// Set the primary variable combination for the 2pnc model
template<class TypeTag>
struct Formulation<TypeTag, TTag::FuelCell>
{ static constexpr auto value = TwoPFormulation::p1s0; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FuelCell>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::H2ON2O2<Scalar>;
};

} // end namespace Dumux::Properties

#endif
