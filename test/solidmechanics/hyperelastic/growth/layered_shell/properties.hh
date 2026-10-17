// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Properties of the growing layered shell test
 */
#ifndef DUMUX_TEST_HYPERELASTIC_GROWTH_LAYERED_SHELL_PROPERTIES_HH
#define DUMUX_TEST_HYPERELASTIC_GROWTH_LAYERED_SHELL_PROPERTIES_HH

#include <dune/alugrid/grid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/solidmechanics/hyperelastic/model.hh>

#include "problem.hh"

#ifndef DIMENSION
#define DIMENSION 2
#endif

namespace Dumux::Properties {

namespace TTag {
struct LayeredShell { using InheritsFrom = std::tuple<Hyperelastic>; };
struct LayeredShellBox { using InheritsFrom = std::tuple<LayeredShell, BoxModel>; };
struct LayeredShellPQ2Hybrid { using InheritsFrom = std::tuple<LayeredShell, PQ2HybridModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::LayeredShell>
{ using type = Dune::ALUGrid<DIMENSION, DIMENSION, Dune::cube, Dune::nonconforming>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::LayeredShell>
{ static constexpr bool value = true; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::LayeredShell>
{
    using type = HyperelasticSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                           GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct Problem<TypeTag, TTag::LayeredShell>
{ using type = LayeredShellProblem<TypeTag>; };

} // end namespace Dumux::Properties

#endif
