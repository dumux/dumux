// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief The properties of the test for the staggered grid Navier-Stokes model with analytical solution (Kovasznay 1948, \cite Kovasznay1948)
 */
#ifndef DUMUX_KOVASZNAY_TEST_PROPERTIES_HH
#define DUMUX_KOVASZNAY_TEST_PROPERTIES_HH

#ifndef UPWINDSCHEMEORDER
#define UPWINDSCHEMEORDER 0
#endif

#include <dune/grid/yaspgrid.hh>

#if HAVE_DUNE_SUBGRID
#include <dune/subgrid/subgrid.hh>
#endif

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct KovasznayTest { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::KovasznayTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::KovasznayTest>
{
    using HostGrid = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >;

#if HAVE_DUNE_SUBGRID
    using type = Dune::SubGrid<HostGrid::dimension, HostGrid>;
#else
    using type = HostGrid;
#endif
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::KovasznayTest> { using type = Dumux::KovasznayTestProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::KovasznayTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::KovasznayTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::KovasznayTest> { static constexpr bool value = true; };

template<class TypeTag>
struct UpwindSchemeOrder<TypeTag, TTag::KovasznayTest> { static constexpr int value = UPWINDSCHEMEORDER; };

} // end namespace Dumux::Properties

#endif
