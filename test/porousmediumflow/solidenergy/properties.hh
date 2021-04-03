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
/**
 * \file
 * \ingroup SolidEnergyTests
 * \brief Properties for the solid energy test
 */
#ifndef DUMUX_TEST_SOLIDENERGY_PROPERTIES_HH
#define DUMUX_TEST_SOLIDENERGY_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/solidenergy/model.hh>

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct SolidEnergyTest { using InheritsFrom = std::tuple<SolidEnergy, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::SolidEnergyTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::SolidEnergyTest> { using type = SolidEnergyProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::SolidEnergyTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = SolidEnergySpatialParams<GridGeometry, Scalar>;
};

} // end namespace Dumux::Properties

#endif
