// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Properties of the dynamic hyperelastic manufactured solution test
 */
#ifndef DUMUX_TEST_HYPERELASTIC_DYNAMIC_MMS_PROPERTIES_HH
#define DUMUX_TEST_HYPERELASTIC_DYNAMIC_MMS_PROPERTIES_HH

#include <dune/alugrid/grid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/solidmechanics/hyperelastic/model.hh>

#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct HyperelasticDynamicMMS { using InheritsFrom = std::tuple<Hyperelastic>; };
struct HyperelasticDynamicMMSPQ2FE { using InheritsFrom = std::tuple<HyperelasticDynamicMMS, PQ2FEModel>; };
struct HyperelasticDynamicMMSPQ2Hybrid { using InheritsFrom = std::tuple<HyperelasticDynamicMMS, PQ2HybridModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::HyperelasticDynamicMMS>
{ using type = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::HyperelasticDynamicMMS>
{ static constexpr bool value = true; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::HyperelasticDynamicMMS>
{
    using type = HyperelasticSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                           GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct Problem<TypeTag, TTag::HyperelasticDynamicMMS>
{ using type = HyperelasticDynamicMMSProblem<TypeTag>; };

} // end namespace Dumux::Properties

#endif
