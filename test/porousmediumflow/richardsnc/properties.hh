// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsNCTests
 * \brief The properties of the water infiltration problem with a
 *        low-permeability lens embedded into a high-permeability domain which
 *        uses the Richards box model.
 */
#ifndef DUMUX_DONEA_TEST_PROPERTIES_HH
#define DUMUX_DONEA_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/richardsnc/model.hh>

#include "spatialparams.hh"
#include "problem.hh"

// Specify the properties for the lens problem
namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RichardsWellTracer { using InheritsFrom = std::tuple<RichardsNC>; };
struct RichardsWellTracerBox { using InheritsFrom = std::tuple<RichardsWellTracer, BoxModel>; };
struct RichardsWellTracerCC { using InheritsFrom = std::tuple<RichardsWellTracer, CCTpfaModel>; };
} // end namespace TTag

// Use 2d YaspGrid
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsWellTracer> { using type = Dune::YaspGrid<2>; };

// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsWellTracer> { using type = RichardsWellTracerProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsWellTracer>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RichardsWellTracerSpatialParams<GridGeometry, Scalar>;
};

// Set the physical problem to be solved
template<class TypeTag>
struct PointSource<TypeTag, TTag::RichardsWellTracer> { using type = SolDependentPointSource<TypeTag>; };

} // end namespace Dumux::Properties

#endif
