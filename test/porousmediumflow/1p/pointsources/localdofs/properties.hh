// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief The properties of the point source tests for control-volume, hybrid and finite element schemes
 */
#ifndef DUMUX_ONEP_POINTSOURCE_LOCALDOFS_PROPERTIES_HH
#define DUMUX_ONEP_POINTSOURCE_LOCALDOFS_PROPERTIES_HH

#include <dune/grid/onedgrid.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1.hh>
#include <dumux/discretization/pq2.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/porousmediumflow/1p/model.hh>

#include "problem.hh"
#include "spatialparams.hh"

#ifndef GRIDTYPE
#define GRIDTYPE Dune::YaspGrid<2>
#endif

namespace Dumux::Properties {

namespace TTag {
struct OnePPointSource { using InheritsFrom = std::tuple<OneP>; };
struct OnePPointSourceBox { using InheritsFrom = std::tuple<OnePPointSource, BoxModel>; };
struct OnePPointSourcePQ2Hybrid { using InheritsFrom = std::tuple<OnePPointSource, PQ2HybridModel>; };
struct OnePPointSourcePQ1FE { using InheritsFrom = std::tuple<OnePPointSource, PQ1FEModel>; };
struct OnePPointSourcePQ2FE { using InheritsFrom = std::tuple<OnePPointSource, PQ2FEModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::OnePPointSource> { using type = GRIDTYPE; };

template<class TypeTag>
struct Problem<TypeTag, TTag::OnePPointSource>
{ using type = OnePPointSourceProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePPointSource>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::Constant<1, Scalar>>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePPointSource>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePPointSourceSpatialParams<GridGeometry, Scalar>;
};

//! Assemble with the variables per local dof
template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::OnePPointSourceBox> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePPointSourcePQ1FE> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePPointSourcePQ2FE> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
