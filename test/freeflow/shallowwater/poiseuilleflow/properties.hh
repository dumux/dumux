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

#ifndef DUMUX_POISEUILLE_FLOW_TEST_PROPERTIES_HH
#define DUMUX_POISEUILLE_FLOW_TEST_PROPERTIES_HH

#include <dumux/freeflow/shallowwater/model.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#ifndef GRIDTYPE
#define GRIDTYPE Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>
#endif

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

namespace TTag {
struct PoiseuilleFlow { using InheritsFrom = std::tuple<ShallowWater, CCTpfaModel>; };
} // namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::PoiseuilleFlow> { using type = GRIDTYPE; };

template<class TypeTag>
struct Problem<TypeTag, TTag::PoiseuilleFlow> { using type = Dumux::PoiseuilleFlowProblem<TypeTag> ; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PoiseuilleFlow>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;

public:
    using type = PoiseuilleFlowSpatialParams<GridGeometry, Scalar, VolumeVariables>;
};

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PoiseuilleFlow> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::PoiseuilleFlow> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
