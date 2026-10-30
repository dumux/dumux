// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief The properties of the single-phase test with two layers of different permeability
 */
#ifndef DUMUX_ONEP_LAYERED_TEST_PROPERTIES_HH
#define DUMUX_ONEP_LAYERED_TEST_PROPERTIES_HH

#include <type_traits>

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/concepts.hh>
#include <dumux/discretization/gridvariables.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>

#include <dumux/flux/cvfe/darcyslaw_.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/1p/variables.hh>
#include <dumux/porousmediumflow/immiscible/localresidual_.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "spatialparams.hh"
#include "problem.hh"

#ifndef GRIDTYPE
#define GRIDTYPE Dune::YaspGrid<2>
#endif

namespace Dumux::Properties {

namespace TTag {
struct OnePLayered { using InheritsFrom = std::tuple<OneP>; };
struct OnePLayeredBox { using InheritsFrom = std::tuple<OnePLayered, BoxModel>; };
struct OnePLayeredPQ1Bubble { using InheritsFrom = std::tuple<OnePLayered, PQ1BubbleModel>; };
struct OnePLayeredPQ1BubbleHybrid { using InheritsFrom = std::tuple<OnePLayered, PQ1BubbleHybridModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::OnePLayered> { using type = GRIDTYPE; };

template<class TypeTag>
struct Problem<TypeTag, TTag::OnePLayered> { using type = OnePLayeredProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePLayered>
{
    using type = OnePLayeredSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                          GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePLayered>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
};

template<class TypeTag>
struct GridVariables<TypeTag, TTag::OnePLayered>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GG::LocalView;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using VariablesTraits = OnePVolumeVariablesTraits<
        GetPropType<TypeTag, Properties::PrimaryVariables>,
        GetPropType<TypeTag, Properties::FluidSystem>,
        GetPropType<TypeTag, Properties::FluidState>,
        GetPropType<TypeTag, Properties::SolidSystem>,
        GetPropType<TypeTag, Properties::SolidState>,
        typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType,
        GetPropType<TypeTag, Properties::ModelTraits>
    >;
    using Variables = Experimental::OnePVariables<VariablesTraits>;
    using IPDataCache = CVFE::LocalBasisInterpolationPointData<GG>;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();

    using HybridCache = Experimental::CVFE::HybridCVFEGridVariablesCache<
        Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>, enableCache>;
    using CVFECache = Experimental::CVFE::CVFEGridVariablesCache<
        Experimental::CVFE::CVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>, enableCache>;
    using Cache = std::conditional_t<Experimental::Concepts::HybridElementDiscretization<ElementDiscretization>,
                                     HybridCache, CVFECache>;
public:
    using type = Experimental::GridVariables<GG, Cache>;
};

template<class TypeTag>
struct AdvectionType<TypeTag, TTag::OnePLayered>
{
    using type = Experimental::CVFEDarcysLawAtIp<GetPropType<TypeTag, Properties::Scalar>,
                                                 GetPropType<TypeTag, Properties::GridGeometry>>;
};

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePLayered>
{ using type = Experimental::ImmiscibleLocalResidual<TypeTag>; };

} // end namespace Dumux::Properties

#endif
