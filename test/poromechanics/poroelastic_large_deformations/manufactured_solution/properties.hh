// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Properties of the plane-strain manufactured solution for the large-deformation
 *        poroelastic model: quadratic displacement (PQ2 hybrid, or PQ2 finite elements with
 *        MOMENTUM_PQ2FE), linear pressures (box).
 */
#ifndef DUMUX_TEST_POROELASTIC_LARGE_DEF_MMS_PROPERTIES_HH
#define DUMUX_TEST_POROELASTIC_LARGE_DEF_MMS_PROPERTIES_HH

#include <tuple>

#include <dune/alugrid/grid.hh>

#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/fem/gridvariablescache.hh>
#include <dumux/discretization/box.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/discretization/gridvariables.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>

#include <dumux/poromechanics/poroelastic_large_deformations/model.hh>
#include <dumux/poromechanics/poroelastic_large_deformations/couplingmanager.hh>

#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct PoroElasticLargeDefMMS { using Grid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>; };
#if MOMENTUM_PQ2FE
struct PoroElasticLargeDefMMSMomentum { using InheritsFrom = std::tuple<PoroElasticLargeDefMMS, PoroElasticLargeDefMomentumModel, PQ2FEModel>; };
#else
struct PoroElasticLargeDefMMSMomentum { using InheritsFrom = std::tuple<PoroElasticLargeDefMMS, PoroElasticLargeDefMomentumModel, PQ2HybridModel>; };
#endif
struct PoroElasticLargeDefMMSSolidPressure { using InheritsFrom = std::tuple<PoroElasticLargeDefMMS, PoroElasticLargeDefSolidPressureModel, BoxModel>; };
struct PoroElasticLargeDefMMSFluidPressure { using InheritsFrom = std::tuple<PoroElasticLargeDefMMS, PoroElasticLargeDefFluidPressureModel, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PoroElasticLargeDefMMS>
{
    using type = PoroElasticLargeDefMMSSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                                     GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::PoroElasticLargeDefMMS> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PoroElasticLargeDefMMS> { static constexpr bool value = true; };

template<class TypeTag>
struct Problem<TypeTag, TTag::PoroElasticLargeDefMMSMomentum>
{ using type = PoroElasticLargeDefMMSMomentumProblem<TypeTag>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::PoroElasticLargeDefMMSSolidPressure>
{ using type = PoroElasticLargeDefMMSSolidPressureProblem<TypeTag>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::PoroElasticLargeDefMMSFluidPressure>
{ using type = PoroElasticLargeDefMMSFluidPressureProblem<TypeTag>; };

template<class TypeTag>
struct GridVariables<TypeTag, TTag::PoroElasticLargeDefMMSMomentum>
{
private:
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
#if MOMENTUM_PQ2FE
    using Traits = Experimental::FEDefaultGridVariablesCacheTraits<
        GetPropType<TypeTag, Properties::Problem>,
        GetPropType<TypeTag, Properties::VolumeVariables>,
        CVFE::LocalBasisInterpolationPointData<GridDiscretization>>;
    using Cache = Experimental::FEGridVariablesCache<Traits, true>;
#else
    using Traits = Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<
        GetPropType<TypeTag, Properties::Problem>,
        GetPropType<TypeTag, Properties::VolumeVariables>,
        CVFE::LocalBasisInterpolationPointData<GridDiscretization>>;
    using Cache = Experimental::CVFE::HybridCVFEGridVariablesCache<Traits, true>;
#endif
public:
    using type = Experimental::GridVariables<GridDiscretization, Cache>;
};

template<class TypeTag>
struct PoroElasticLargeDefMMSBoxGridVariables
{
private:
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using Traits = Experimental::CVFE::CVFEDefaultGridVariablesCacheTraits<
        GetPropType<TypeTag, Properties::Problem>,
        GetPropType<TypeTag, Properties::VolumeVariables>,
        CVFE::LocalBasisInterpolationPointData<GridDiscretization>>;
    using Cache = Experimental::CVFE::CVFEGridVariablesCache<Traits, true>;
public:
    using type = Experimental::GridVariables<GridDiscretization, Cache>;
};

template<class TypeTag>
struct GridVariables<TypeTag, TTag::PoroElasticLargeDefMMSSolidPressure>
{ using type = typename PoroElasticLargeDefMMSBoxGridVariables<TypeTag>::type; };

template<class TypeTag>
struct GridVariables<TypeTag, TTag::PoroElasticLargeDefMMSFluidPressure>
{ using type = typename PoroElasticLargeDefMMSBoxGridVariables<TypeTag>::type; };

using PoroElasticLargeDefMMSTraits = MultiDomainTraits<
    TTag::PoroElasticLargeDefMMSMomentum,
    TTag::PoroElasticLargeDefMMSSolidPressure,
    TTag::PoroElasticLargeDefMMSFluidPressure>;

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PoroElasticLargeDefMMS>
{ using type = PoroElasticLargeDefCouplingManager<PoroElasticLargeDefMMSTraits>; };

} // end namespace Dumux::Properties

#endif
