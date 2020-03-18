// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief The properties for the shallow water bowl test
 */
#ifndef DUMUX_BOWL_TEST_PROPERTIES_HH
#define DUMUX_BOWL_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/common/properties.hh>
#include <dumux/freeflow/shallowwater/model.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct Bowl { using InheritsFrom = std::tuple<ShallowWater, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::Bowl>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Bowl>
{ using type = Dumux::BowlProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Bowl>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
public:
    using type = BowlSpatialParams<GridGeometry, Scalar, VolumeVariables>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Bowl>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Bowl>
{ static constexpr bool value = false; };

} // end namespace Dumux::Properties

#endif
