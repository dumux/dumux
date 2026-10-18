// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Properties of the ring collapse test
 */
#ifndef DUMUX_TEST_CONTACT_RING_PROPERTIES_HH
#define DUMUX_TEST_CONTACT_RING_PROPERTIES_HH

#include <dune/alugrid/grid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/box.hh>
#include <dumux/solidmechanics/hyperelastic/model.hh>
#include <dumux/solidmechanics/hyperelastic/spatialparams.hh>

#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct RingCollapse { using InheritsFrom = std::tuple<Hyperelastic, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::RingCollapse>
{ using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::RingCollapse>
{ static constexpr bool value = true; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RingCollapse>
{
    using type = HyperelasticSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                           GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct Problem<TypeTag, TTag::RingCollapse>
{ using type = RingCollapseProblem<TypeTag>; };

} // end namespace Dumux::Properties

#endif
