// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// SPP-1748 incompressible block mixed formulations:
// - 3D P1BP1 (mini/mixed, PQ1Bubble displacement + P1 pressure)
// - 3D P2P1 (Taylor-Hood, PQ2 displacement + P1 pressure)
#ifndef DUMUX_HYPERELASTIC_INCOMPRESSIBLE_BLOCK_PROPERTIES_HH
#define DUMUX_HYPERELASTIC_INCOMPRESSIBLE_BLOCK_PROPERTIES_HH

#include <type_traits>

#include <dune/alugrid/grid.hh>

#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/box.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/discretization/gridvariables.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>

#include <dumux/solidmechanics/hyperelastic_voliso/model.hh>
#include <dumux/solidmechanics/hyperelastic_voliso/couplingmanager.hh>
#include "../../hyperelastic/incompressible_block/spatialparams.hh"
#include "problem.hh"

#ifndef ALU_ELEMENT_TYPE
#define ALU_ELEMENT_TYPE Dune::simplex
#endif

namespace Dumux::Properties {

namespace TTag {

struct IncompressibleBlockCommon {};

struct IncompressibleBlockMomentumP1B
{
    template<class TypeTag>
    using Problem = IncompressibleBlockMomentumProblem<TypeTag>;

    using Grid = Dune::ALUGrid<3, 3, ALU_ELEMENT_TYPE, Dune::nonconforming>;

    using InheritsFrom = std::tuple<
        HyperelasticVolIsoMomentumModel,
        PQ1BubbleModel,
        IncompressibleBlockCommon>;

    using EnableGridGeometryCache = std::true_type;
    using EnableGridVolumeVariablesCache = std::false_type;
};

struct IncompressibleBlockMomentumPQ2
{
    template<class TypeTag>
    using Problem = IncompressibleBlockMomentumProblem<TypeTag>;

    using Grid = Dune::ALUGrid<3, 3, ALU_ELEMENT_TYPE, Dune::nonconforming>;

    using InheritsFrom = std::tuple<
        HyperelasticVolIsoMomentumModel,
        PQ2HybridModel,
        IncompressibleBlockCommon>;

    using EnableGridGeometryCache = std::true_type;
    using EnableGridVolumeVariablesCache = std::false_type;
};

struct IncompressibleBlockMomentumPQ2Fe
{
    template<class TypeTag>
    using Problem = IncompressibleBlockMomentumProblem<TypeTag>;

    using Grid = Dune::ALUGrid<3, 3, ALU_ELEMENT_TYPE, Dune::nonconforming>;

    using InheritsFrom = std::tuple<
        HyperelasticVolIsoMomentumModel,
        PQ2FEModel,
        IncompressibleBlockCommon>;

    using EnableGridGeometryCache = std::true_type;
    using EnableGridVolumeVariablesCache = std::false_type;
};

struct IncompressibleBlockPressure
{
    template<class TypeTag>
    using Problem = IncompressibleBlockPressureProblem<TypeTag>;

    using Grid = Dune::ALUGrid<3, 3, ALU_ELEMENT_TYPE, Dune::nonconforming>;

    using InheritsFrom = std::tuple<
        HyperelasticVolIsoPressureModel,
        BoxModel,
        IncompressibleBlockCommon>;

    using EnableGridGeometryCache = std::true_type;
    using EnableGridVolumeVariablesCache = std::false_type;
};

struct IncompressibleBlockPressureFe
{
    template<class TypeTag>
    using Problem = IncompressibleBlockPressureProblem<TypeTag>;

    using Grid = Dune::ALUGrid<3, 3, ALU_ELEMENT_TYPE, Dune::nonconforming>;

    using InheritsFrom = std::tuple<
        HyperelasticVolIsoPressureModel,
        PQ1FEModel,
        IncompressibleBlockCommon>;

    using EnableGridGeometryCache = std::true_type;
    using EnableGridVolumeVariablesCache = std::false_type;
};

} // end namespace TTag

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::IncompressibleBlockCommon>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = IncompressibleBlockSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct GridVariables<TypeTag, TTag::IncompressibleBlockMomentumP1B>
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
struct GridVariables<TypeTag, TTag::IncompressibleBlockMomentumPQ2>
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
struct GridVariables<TypeTag, TTag::IncompressibleBlockPressure>
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

#ifndef MOMENTUM_DISCRETIZATION
#define MOMENTUM_DISCRETIZATION IncompressibleBlockMomentumP1B
#endif

#ifndef PRESSURE_DISCRETIZATION
#define PRESSURE_DISCRETIZATION IncompressibleBlockPressure
#endif

template<class TypeTagA, class TypeTagB>
using TheCouplingManager = HyperelasticVolIsoCouplingManager<MultiDomainTraits<TypeTagA, TypeTagB>>;

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::MOMENTUM_DISCRETIZATION>
{ using type = TheCouplingManager<TTag::MOMENTUM_DISCRETIZATION, TTag::PRESSURE_DISCRETIZATION>; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PRESSURE_DISCRETIZATION>
{ using type = TheCouplingManager<TTag::MOMENTUM_DISCRETIZATION, TTag::PRESSURE_DISCRETIZATION>; };

} // end namespace Dumux::Properties
#endif