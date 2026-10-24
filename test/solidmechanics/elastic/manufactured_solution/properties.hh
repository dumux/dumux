// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsTests
 * \brief Properties of the linear elastic manufactured solution, assembled with local degrees of freedom
 */
#ifndef DUMUX_ELASTIC_PROPERTIES_HH
#define DUMUX_ELASTIC_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/solidmechanics/elastic/model.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct TestElastic { using InheritsFrom = std::tuple<Elastic, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::TestElastic> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::TestElastic> { using type = Dumux::ElasticProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TestElastic>
{
    using type = ElasticManufacturedSolutionSpatialParams<GetPropType<TypeTag, Properties::Scalar>,
                                                          GetPropType<TypeTag, Properties::GridGeometry>>;
};

//! Assemble with the variables per local dof
template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::TestElastic> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
