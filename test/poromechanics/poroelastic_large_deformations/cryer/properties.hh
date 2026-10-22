// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Properties of the Cryer sphere consolidation benchmark: quadratic displacement
 *        (PQ2 hybrid), linear pressures (box).
 */
#ifndef DUMUX_TEST_POROELASTIC_LARGE_DEF_CRYER_PROPERTIES_HH
#define DUMUX_TEST_POROELASTIC_LARGE_DEF_CRYER_PROPERTIES_HH

#include <tuple>

#include <dune/alugrid/grid.hh>

#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/box.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/poromechanics/poroelastic_large_deformations/model.hh>
#include <dumux/poromechanics/poroelastic_large_deformations/couplingmanager.hh>

#include "spatialparams.hh"
#include "problem_momentum.hh"
#include "problem_solid_pressure.hh"
#include "problem_fluid_pressure.hh"

namespace Dumux::Properties {

using CryerGrid = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;

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
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = CryerSpatialParams<GridDiscretization, Scalar>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::CryerMomentum>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::CryerMomentum>
{ static constexpr bool value = true; };

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
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = CryerSpatialParams<GridDiscretization, Scalar>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::CryerSolidPressure>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::CryerSolidPressure>
{ static constexpr bool value = true; };

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
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = CryerSpatialParams<GridDiscretization, Scalar>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::CryerFluidPressure>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::CryerFluidPressure>
{ static constexpr bool value = true; };

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
