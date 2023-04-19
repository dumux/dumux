// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MPNCTests
 * \brief The properties of the problem where air is injected in a unsaturated porous medium.
 */
#ifndef DUMUX_MPNC_TWOPTWOC_COMPARISON_OBSTACLE_PROPERTIES_HH
#define DUMUX_MPNC_TWOPTWOC_COMPARISON_OBSTACLE_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/mpnc/model.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include <test/porousmediumflow/2p2c/mpnccomparison/iofields.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct MPNCComparison { using InheritsFrom = std::tuple<MPNC>; };
struct MPNCComparisonBox { using InheritsFrom = std::tuple<MPNCComparison, BoxModel>; };
struct MPNCComparisonCC { using InheritsFrom = std::tuple<MPNCComparison, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::MPNCComparison> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::MPNCComparison> { using type = MPNCComparisonProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::MPNCComparison>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = MPNCComparisonSpatialParams<GridGeometry, Scalar>;
};

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::MPNCComparison>
{
    using type = FluidSystems::H2ON2<GetPropType<TypeTag, Properties::Scalar>,
                                     FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>;
};

template<class TypeTag>
struct UseMoles<TypeTag, TTag::MPNCComparison> { static constexpr bool value = true; };
template<class TypeTag>
struct IOFields<TypeTag, TTag::MPNCComparison> { using type = TwoPTwoCMPNCIOFields; };

} // end namespace Dumux::Properties

#endif
