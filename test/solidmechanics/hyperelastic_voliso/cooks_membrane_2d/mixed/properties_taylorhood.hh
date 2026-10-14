// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_COOKS_MEMBRANE_TAYLORHOOD_PROPERTIES_HH
#define DUMUX_COOKS_MEMBRANE_TAYLORHOOD_PROPERTIES_HH

#include <type_traits>

#include <dune/alugrid/grid.hh>

#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/box.hh>
#include <dumux/multidomain/traits.hh>

// New variables concept (hybrid for PQ2 momentum, standard for Box pressure)
#include <dumux/discretization/gridvariables.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>

#include <dumux/solidmechanics/hyperelastic_voliso/model.hh>
#include <dumux/solidmechanics/hyperelastic_voliso/couplingmanager.hh>
#include "../../../hyperelastic/cooks_membrane_2d/spatialparams.hh"
#include "problem_momentum_mixed.hh"
#include "problem_pressure_mixed.hh"

namespace Dumux::Properties {

// ---- Momentum type tag (PQ2 hybrid displacement — Taylor-Hood) ----
namespace TTag {
struct CooksMembraneTHMomentum
{
    using InheritsFrom = std::tuple<HyperelasticVolIsoMomentumModel, PQ2HybridModel>;
    using Grid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>;
};
} // end TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::CooksMembraneTHMomentum>
{ using type = CooksMembraneMomentumMixedProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CooksMembraneTHMomentum>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = CooksMembraneHyperelasticSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::CooksMembraneTHMomentum>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::CooksMembraneTHMomentum>
{ static constexpr bool value = true; };

//! Taylor-Hood uses HybridCVFEGridVariablesCache (stores element-level QP data for edge DOFs).
template<class TypeTag>
struct GridVariables<TypeTag, TTag::CooksMembraneTHMomentum>
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

// ---- Pressure type tag (Box/P1 bulk pressure — same as MINI) ----
namespace TTag {
struct CooksMembraneTHPressure
{
    using InheritsFrom = std::tuple<HyperelasticVolIsoPressureModel, BoxModel>;
    using Grid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>;
};
} // end TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::CooksMembraneTHPressure>
{ using type = CooksMembranePressureMixedProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CooksMembraneTHPressure>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = CooksMembraneHyperelasticSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::CooksMembraneTHPressure>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::CooksMembraneTHPressure>
{ static constexpr bool value = true; };

template<class TypeTag>
struct GridVariables<TypeTag, TTag::CooksMembraneTHPressure>
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
using CooksMembraneTHMDTraits = MultiDomainTraits<
    TTag::CooksMembraneTHMomentum,
    TTag::CooksMembraneTHPressure>;

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::CooksMembraneTHMomentum>
{ using type = HyperelasticVolIsoCouplingManager<CooksMembraneTHMDTraits>; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::CooksMembraneTHPressure>
{ using type = HyperelasticVolIsoCouplingManager<CooksMembraneTHMDTraits>; };

} // end namespace Dumux::Properties
#endif
