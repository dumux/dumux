// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief The properties for the convergence test with analytic solution
 */
#ifndef DUMUX_CONVERGENCE_TEST_ONEP_PROPERTIES_HH
#define DUMUX_CONVERGENCE_TEST_ONEP_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include "spatialparams.hh"
#include "problem.hh"

#ifndef GRIDTYPE
#define GRIDTYPE Dune::YaspGrid<2>
#endif

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct OnePConvergence { using InheritsFrom = std::tuple<OneP>; };
struct OnePConvergenceTpfa { using InheritsFrom = std::tuple<OnePConvergence, CCTpfaModel>; };
struct OnePConvergenceMpfa { using InheritsFrom = std::tuple<OnePConvergence, CCMpfaModel>; };
struct OnePConvergenceBox { using InheritsFrom = std::tuple<OnePConvergence, BoxModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePConvergence> { using type = Dumux::ConvergenceProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePConvergence>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::Constant<1, Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePConvergence> { using type = GRIDTYPE; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePConvergence>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = ConvergenceTestSpatialParams<GridGeometry, Scalar>;
};

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePConvergence> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePConvergence> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePConvergence> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
