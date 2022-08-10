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
 * \ingroup NavierStokesTests
 * \brief the properties of the freeflow problem for pipe flow
 * Simulation of a radially-symmetric pipe flow with circular cross-section
 */
#ifndef DUMUX_TEST_FREEFLOW_PIPE_PROPERTIES_HH
#define DUMUX_TEST_FREEFLOW_PIPE_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/freeflow/navierstokes/momentum/diamond/model.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>

#include <dumux/discretization/fcdiamond.hh>
#include <dumux/discretization/cvfe.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct PipeFlow { };
struct PipeFlowDiamond { using InheritsFrom = std::tuple<NavierStokesMomentumDiamond, FaceCenteredDiamondModel, PipeFlow>; };
struct PipeFlowCvfe { using InheritsFrom = std::tuple<NavierStokesMomentumCvfe, CvfeModel, PipeFlow>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PipeFlow>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::Constant<1, Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PipeFlow>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PipeFlow>
{ using type = FreeFlowPipeProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::PipeFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::PipeFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PipeFlow> { static constexpr bool value = true; };

// rotation-symmetric grid geometry forming a cylinder channel
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PipeFlowDiamond>
{
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;

    struct GGTraits : public FaceCenteredDiamondDefaultGridGeometryTraits<GridView>
    { using Extrusion = RotationalExtrusion<0>; };

    using type = FaceCenteredDiamondFVGridGeometry<GridView, enableCache, GGTraits>;
};

// rotation-symmetric grid geometry forming a cylinder channel
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PipeFlowCvfe>
{
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    struct GGTraits : public CvfeDefaultGridGeometryTraits<GridView>
    { using Extrusion = RotationalExtrusion<0>; };

    using type = CvfeFVGridGeometry<Scalar, GridView, enableCache, GGTraits>;
};

} // end namespace Dumux::Properties

#endif
