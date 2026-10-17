// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Properties of the growth-induced wrinkling test
 */
#ifndef DUMUX_TEST_HYPERELASTIC_GROWTH_WRINKLING_PROPERTIES_HH
#define DUMUX_TEST_HYPERELASTIC_GROWTH_WRINKLING_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/solidmechanics/hyperelastic/model.hh>

#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct FilmSubstrate { using InheritsFrom = std::tuple<Hyperelastic>; };
struct FilmSubstrateBox { using InheritsFrom = std::tuple<FilmSubstrate, BoxModel>; };
struct FilmSubstratePQ2Hybrid { using InheritsFrom = std::tuple<FilmSubstrate, PQ2HybridModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::FilmSubstrate>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<double, 2>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::FilmSubstrate>
{ static constexpr bool value = true; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::FilmSubstrate>
{
    using type = HyperelasticSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                           GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct Problem<TypeTag, TTag::FilmSubstrate>
{ using type = FilmSubstrateProblem<TypeTag>; };

} // end namespace Dumux::Properties

#endif
