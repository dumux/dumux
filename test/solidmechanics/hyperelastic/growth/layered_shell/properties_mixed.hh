// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Properties of the growing layered shell test in the mixed displacement-pressure formulation
 */
#ifndef DUMUX_TEST_HYPERELASTIC_GROWTH_LAYERED_SHELL_PROPERTIES_MIXED_HH
#define DUMUX_TEST_HYPERELASTIC_GROWTH_LAYERED_SHELL_PROPERTIES_MIXED_HH

#include <type_traits>

#include <dune/alugrid/grid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/gridvariables.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/solidmechanics/hyperelastic/spatialparams.hh>
#include <dumux/solidmechanics/hyperelastic_voliso/couplingmanager.hh>
#include <dumux/solidmechanics/hyperelastic_voliso/model.hh>

#include "problem_mixed.hh"

#ifndef DIMENSION
#define DIMENSION 2
#endif

namespace Dumux::Properties {

namespace TTag {
struct LayeredShellMixed
{
    using Grid = Dune::ALUGrid<DIMENSION, DIMENSION, Dune::cube, Dune::nonconforming>;
    using EnableGridGeometryCache = std::true_type;
    using EnableGridVolumeVariablesCache = std::false_type;
};

struct LayeredShellMomentumP1B
{
    using InheritsFrom = std::tuple<HyperelasticVolIsoMomentumModel, PQ1BubbleModel, LayeredShellMixed>;
    template<class TypeTag>
    using Problem = LayeredShellMixedMomentumProblem<TypeTag>;
};

struct LayeredShellMomentumQ2
{
    using InheritsFrom = std::tuple<HyperelasticVolIsoMomentumModel, PQ2HybridModel, LayeredShellMixed>;
    template<class TypeTag>
    using Problem = LayeredShellMixedMomentumProblem<TypeTag>;
};

struct LayeredShellPressure
{
    using InheritsFrom = std::tuple<HyperelasticVolIsoPressureModel, BoxModel, LayeredShellMixed>;
    template<class TypeTag>
    using Problem = LayeredShellMixedPressureProblem<TypeTag>;
};
} // end namespace TTag

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::LayeredShellMixed>
{
    using type = HyperelasticSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                           GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct GridVariables<TypeTag, TTag::LayeredShellMomentumP1B>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Variables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::CVFE::CVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::CVFE::CVFEGridVariablesCache<Traits, enableCache>;
public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

template<class TypeTag>
struct GridVariables<TypeTag, TTag::LayeredShellMomentumQ2>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Variables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::CVFE::HybridCVFEGridVariablesCache<Traits, enableCache>;
public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

template<class TypeTag>
struct GridVariables<TypeTag, TTag::LayeredShellPressure>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Variables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::CVFE::CVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::CVFE::CVFEGridVariablesCache<Traits, enableCache>;
public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

#ifndef MOMENTUM_TYPETAG
#define MOMENTUM_TYPETAG LayeredShellMomentumP1B
#endif

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::MOMENTUM_TYPETAG>
{ using type = HyperelasticVolIsoCouplingManager<MultiDomainTraits<TTag::MOMENTUM_TYPETAG, TTag::LayeredShellPressure>>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::LayeredShellPressure>
{ using type = HyperelasticVolIsoCouplingManager<MultiDomainTraits<TTag::MOMENTUM_TYPETAG, TTag::LayeredShellPressure>>; };

} // end namespace Dumux::Properties

#endif
