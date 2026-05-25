// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_COOKS_MEMBRANE_MIXED_PROPERTIES_HH
#define DUMUX_COOKS_MEMBRANE_MIXED_PROPERTIES_HH

#include <type_traits>

#include <dune/alugrid/grid.hh>

#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/box.hh>
#include <dumux/multidomain/traits.hh>

// New variables concept
#include <dumux/discretization/gridvariables.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>

#include <dumux/solidmechanics/hyperelastic_voliso/model.hh>
#include <dumux/solidmechanics/hyperelastic_voliso/couplingmanager.hh>
#include "../../../hyperelastic/cooks_membrane_2d/spatialparams.hh"
#include "problem_momentum_mixed.hh"
#include "problem_pressure_mixed.hh"

namespace Dumux::Properties {

// ---- Momentum type tag (PQ1Bubble displacement) ----
namespace TTag {
struct CooksMembraneMixedMomentum
{
    using InheritsFrom = std::tuple<HyperelasticVolIsoMomentumModel, PQ1BubbleModel>;
    using Grid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>;
};
} // end TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::CooksMembraneMixedMomentum>
{ using type = CooksMembraneMomentumMixedProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CooksMembraneMixedMomentum>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = CooksMembraneHyperelasticSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::CooksMembraneMixedMomentum>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::CooksMembraneMixedMomentum>
{ static constexpr bool value = true; };

//! Use Experimental::GridVariables for momentum so that ElementVariables
//! covers all local DOFs (vertices + bubble) via CVFEGridVariablesCache.
template<class TypeTag>
struct GridVariables<TypeTag, TTag::CooksMembraneMixedMomentum>
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

// ---- Pressure type tag (Box/P1 bulk pressure) ----
namespace TTag {
struct CooksMembraneMixedPressure
{
    using InheritsFrom = std::tuple<HyperelasticVolIsoPressureModel, BoxModel>;
    using Grid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>;
};
} // end TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::CooksMembraneMixedPressure>
{ using type = CooksMembranePressureMixedProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CooksMembraneMixedPressure>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = CooksMembraneHyperelasticSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::CooksMembraneMixedPressure>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::CooksMembraneMixedPressure>
{ static constexpr bool value = true; };

//! Use Experimental::GridVariables for pressure as well.
template<class TypeTag>
struct GridVariables<TypeTag, TTag::CooksMembraneMixedPressure>
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

// ---- Coupling manager ----
using CooksMembraneMixedMDTraits = MultiDomainTraits<
    TTag::CooksMembraneMixedMomentum,
    TTag::CooksMembraneMixedPressure>;

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::CooksMembraneMixedMomentum>
{ using type = HyperelasticVolIsoCouplingManager<CooksMembraneMixedMDTraits>; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::CooksMembraneMixedPressure>
{ using type = HyperelasticVolIsoCouplingManager<CooksMembraneMixedMDTraits>; };

} // end namespace Dumux::Properties
#endif
