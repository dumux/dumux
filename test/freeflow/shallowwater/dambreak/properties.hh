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
 * \ingroup ShallowWaterTests
 * \brief The properties of the test for the Shallow water model (wet dam break).
 */
#ifndef DUMUX_DAM_BREAK_TEST_PROPERTIES_HH
#define DUMUX_DAM_BREAK_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/freeflow/shallowwater/model.hh>

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DamBreakWet { using InheritsFrom = std::tuple<ShallowWater, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::DamBreakWet>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DamBreakWet>
{ using type = Dumux::DamBreakProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DamBreakWet>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = DamBreakSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DamBreakWet>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::DamBreakWet>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::DamBreakWet>
{ static constexpr bool value = false; };

} // end namespace Dumux::Properties

#endif
