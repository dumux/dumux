// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief A test problem for the one-phase model:
 * Water is injected in one single point in the middle of the domain.
 */
#ifndef DUMUX_1P_SINGULARITY_PROBLEM_PROPERTIES_HH
#define DUMUX_1P_SINGULARITY_PROBLEM_PROPERTIES_HH

#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/fvspatialparams1pconstant.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#ifndef GRIDTYPE // default to yasp grid if not provided by CMake
#define GRIDTYPE Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >
#endif

#include "problem.hh"
namespace Dumux::Properties {
// Create new type tags
namespace TTag {
struct OnePSingularity { using InheritsFrom = std::tuple<OneP>; };
struct OnePSingularityBox { using InheritsFrom = std::tuple<OnePSingularity, BoxModel>; };
struct OnePSingularityCCTpfa { using InheritsFrom = std::tuple<OnePSingularity, CCTpfaModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePSingularity>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePSingularity> { using type = GRIDTYPE; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePSingularity> { using type = OnePSingularityProblem<TypeTag> ; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePSingularity>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FVPorousMediumFlowSpatialParamsOnePConstant<GridGeometry, Scalar>;
};

} // end namespace Dumux

#endif
