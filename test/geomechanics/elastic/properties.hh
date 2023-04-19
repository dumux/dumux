// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsTests
 * \brief The properties of a test problem for the linear elastic model.
 */
#ifndef DUMUX_ELASTIC_PROPERTIES_HH
#define DUMUX_ELASTIC_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/box.hh>
#include <dumux/geomechanics/elastic/model.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tag
namespace TTag {
struct TestElastic { using InheritsFrom = std::tuple<Elastic, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TestElastic> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TestElastic> { using type = Dumux::ElasticProblem<TypeTag>; };

// The spatial parameters property
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TestElastic>
{
    using type = ElasticSpatialParams< GetPropType<TypeTag, Properties::Scalar>,
                                       GetPropType<TypeTag, Properties::GridGeometry> >;
};

} // end namespace Dumux::Properties

#endif
