// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Properties of the growth bending strip test
 */
#ifndef DUMUX_TEST_HYPERELASTIC_GROWTH_BENDING_STRIP_PROPERTIES_HH
#define DUMUX_TEST_HYPERELASTIC_GROWTH_BENDING_STRIP_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

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
struct BendingStrip { using InheritsFrom = std::tuple<Hyperelastic>; };
struct BendingStripBox { using InheritsFrom = std::tuple<BendingStrip, BoxModel>; };
struct BendingStripPQ2Hybrid { using InheritsFrom = std::tuple<BendingStrip, PQ2HybridModel>; };
struct BendingStripPQ2FE { using InheritsFrom = std::tuple<BendingStrip, PQ2FEModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::BendingStrip>
{ using type = Dune::YaspGrid<DIMENSION>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::BendingStrip>
{ static constexpr bool value = true; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BendingStrip>
{
    using type = HyperelasticSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                           GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct Problem<TypeTag, TTag::BendingStrip>
{ using type = BendingStripProblem<TypeTag>; };

} // end namespace Dumux::Properties

#endif
