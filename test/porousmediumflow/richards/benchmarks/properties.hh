// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsTests
 * \brief Common properties for the Richards benchmarks
 */
#ifndef DUMUX_RICHARDS_BENCHMARKS_PROPERTIES_HH
#define DUMUX_RICHARDS_BENCHMARKS_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/richards/model.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RichardsBenchmark { using InheritsFrom = std::tuple<Richards>; };
struct RichardsBenchmarkCCTpfa { using InheritsFrom = std::tuple<RichardsBenchmark, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsBenchmark>
{ using type = Dune::YaspGrid<1, Dune::TensorProductCoordinates<double, 1>>; };

// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsBenchmark>
{ using type = RichardsBenchmarkProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsBenchmark>
{
    using type = RichardsBenchmarkSpatialParams<
        GetPropType<TypeTag, Properties::GridGeometry>, GetPropType<TypeTag, Properties::Scalar>
    >;
};

} // end namespace Dumux::Properties

#endif
