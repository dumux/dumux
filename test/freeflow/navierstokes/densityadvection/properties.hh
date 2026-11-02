// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief The properties of the test for the advection of a density interface in the CVFE momentum model.
 */
#ifndef DUMUX_DENSITY_ADVECTION_TEST_PROPERTIES_HH
#define DUMUX_DENSITY_ADVECTION_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>

#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq2.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct DensityAdvectionTest {};
struct DensityAdvectionTestPQ1BubbleHybrid { using InheritsFrom = std::tuple<DensityAdvectionTest, NavierStokesMomentumCVFE, PQ1BubbleHybridModel>; };
struct DensityAdvectionTestPQ2Hybrid { using InheritsFrom = std::tuple<DensityAdvectionTest, NavierStokesMomentumCVFE, PQ2HybridModel>; };
} // end namespace TTag

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DensityAdvectionTest>
{
    using type = FluidSystems::OnePLiquid<GetPropType<TypeTag, Properties::Scalar>,
                                          Components::Constant<1, GetPropType<TypeTag, Properties::Scalar>>>;
};

template<class TypeTag>
struct Grid<TypeTag, TTag::DensityAdvectionTest>
{ using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::DensityAdvectionTest>
{ using type = DensityAdvectionTestProblem<TypeTag, Dumux::CVFENavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DensityAdvectionTest>
{
    struct EmptyCouplingManager {};
    using type = EmptyCouplingManager;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DensityAdvectionTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::DensityAdvectionTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::DensityAdvectionTest> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
