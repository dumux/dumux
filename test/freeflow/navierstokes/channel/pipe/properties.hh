// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief the properties of the freeflow problem for pipe flow
 * Simulation of a radially-symmetric pipe flow with circular cross-section
 */
#ifndef DUMUX_TEST_FREEFLOW_PIPE_PROPERTIES_HH
#define DUMUX_TEST_FREEFLOW_PIPE_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct PipeFlow {};
struct PipeFlowMomentum { using InheritsFrom = std::tuple<PipeFlow, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct PipeFlowMass { using InheritsFrom = std::tuple<PipeFlow, NavierStokesMassOneP, CCTpfaModel>; };
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
struct Problem<TypeTag, TTag::PipeFlowMomentum>
{ using type = FreeFlowPipeProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::PipeFlowMass>
{ using type = FreeFlowPipeProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::PipeFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::PipeFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PipeFlow> { static constexpr bool value = true; };

// rotation-symmetric grid geometry forming a cylinder channel
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PipeFlowMomentum>
{
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;

    struct GGTraits : public FaceCenteredStaggeredDefaultGridGeometryTraits<GridView>
    { using Extrusion = RotationalExtrusion<0>; };

    using type = FaceCenteredStaggeredFVGridGeometry<GridView, enableCache, GGTraits>;
};

// rotation-symmetric grid geometry forming a cylinder channel
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PipeFlowMass>
{
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;

    struct GGTraits : public CCTpfaDefaultGridGeometryTraits<GridView>
    { using Extrusion = RotationalExtrusion<0>; };

    using type = CCTpfaFVGridGeometry<GridView, enableCache, GGTraits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PipeFlow>
{
private:
    using Traits = MultiDomainTraits<TTag::PipeFlowMomentum, TTag::PipeFlowMass>;
public:
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
