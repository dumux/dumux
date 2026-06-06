// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_CRYER_PROPERTIES_HH
#define DUMUX_CRYER_PROPERTIES_HH

#include <dune/alugrid/grid.hh>

#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/box.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/discretization/gridvariables.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>

#include <dumux/poromechanics/poroelastic_large_deformations/model.hh>
#include <dumux/poromechanics/poroelastic_large_deformations/couplingmanager.hh>

#include "spatialparams.hh"
#include "problem_momentum.hh"
#include "problem_solid_pressure.hh"
#include "problem_fluid_pressure.hh"

namespace Dumux::Properties {

// ---------------------------------------------------------------------------
// Shared grid type: 3D simplicial ALUGrid
// ---------------------------------------------------------------------------
using CryerGrid = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;

// ---------------------------------------------------------------------------
// Helper: set up Experimental::GridVariables for a CVFE subdomain
// ---------------------------------------------------------------------------
template<class TypeTag>
struct CVFEGridVariablesSetting
{
private:
    using GG       = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem  = GetPropType<TypeTag, Properties::Problem>;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Variables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::CVFE::CVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::CVFE::CVFEGridVariablesCache<Traits, enableCache>;
public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

// ---------------------------------------------------------------------------
// Helper: set up Experimental::GridVariables for a hybrid CVFE subdomain
// (PQ2 displacement: vertex CV dofs + edge-midpoint non-CV/FEM dofs)
// ---------------------------------------------------------------------------
template<class TypeTag>
struct HybridCVFEGridVariablesSetting
{
private:
    using GG       = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem  = GetPropType<TypeTag, Properties::Problem>;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Variables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::CVFE::HybridCVFEGridVariablesCache<Traits, enableCache>;
public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

// ---------------------------------------------------------------------------
// Momentum type tag  (PQ1Bubble displacement + experimental GridVariables)
// ---------------------------------------------------------------------------
namespace TTag {
struct CryerMomentum
{
    using InheritsFrom = std::tuple<PoroElasticLargeDefMomentumModel, PQ2HybridModel>;
    using Grid = CryerGrid;
};
} // end TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::CryerMomentum>
{ using type = CryerMomentumProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CryerMomentum>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = CryerSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::CryerMomentum>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::CryerMomentum>
{ static constexpr bool value = true; };

template<class TypeTag>
struct GridVariables<TypeTag, TTag::CryerMomentum>
{ using type = typename HybridCVFEGridVariablesSetting<TypeTag>::type; };

// ---------------------------------------------------------------------------
// Solid pressure type tag  (Box P1)
// ---------------------------------------------------------------------------
namespace TTag {
struct CryerSolidPressure
{
    using InheritsFrom = std::tuple<PoroElasticLargeDefSolidPressureModel, BoxModel>;
    using Grid = CryerGrid;
};
} // end TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::CryerSolidPressure>
{ using type = CryerSolidPressureProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CryerSolidPressure>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = CryerSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::CryerSolidPressure>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::CryerSolidPressure>
{ static constexpr bool value = true; };

template<class TypeTag>
struct GridVariables<TypeTag, TTag::CryerSolidPressure>
{ using type = typename CVFEGridVariablesSetting<TypeTag>::type; };

// ---------------------------------------------------------------------------
// Fluid pressure type tag  (Box P1)
// ---------------------------------------------------------------------------
namespace TTag {
struct CryerFluidPressure
{
    using InheritsFrom = std::tuple<PoroElasticLargeDefFluidPressureModel, BoxModel>;
    using Grid = CryerGrid;
};
} // end TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::CryerFluidPressure>
{ using type = CryerFluidPressureProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CryerFluidPressure>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = CryerSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::CryerFluidPressure>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::CryerFluidPressure>
{ static constexpr bool value = true; };

template<class TypeTag>
struct GridVariables<TypeTag, TTag::CryerFluidPressure>
{ using type = typename CVFEGridVariablesSetting<TypeTag>::type; };

// ---------------------------------------------------------------------------
// MultiDomain traits + coupling manager
// ---------------------------------------------------------------------------
using CryerMDTraits = MultiDomainTraits<
    TTag::CryerMomentum,
    TTag::CryerSolidPressure,
    TTag::CryerFluidPressure>;

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::CryerMomentum>
{ using type = PoroElasticLargeDefCouplingManager<CryerMDTraits>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::CryerSolidPressure>
{ using type = PoroElasticLargeDefCouplingManager<CryerMDTraits>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::CryerFluidPressure>
{ using type = PoroElasticLargeDefCouplingManager<CryerMDTraits>; };

} // end namespace Dumux::Properties

#endif
