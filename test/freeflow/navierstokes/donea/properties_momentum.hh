// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test for the staggered grid (Navier-)Stokes model with analytical solution (Donea 2003, \cite Donea2003).
 */
#ifndef DUMUX_DONEA_TEST_PROPERTIES_MOMENTUM_HH
#define DUMUX_DONEA_TEST_PROPERTIES_MOMENTUM_HH

#ifndef ENABLECACHING
#define ENABLECACHING true
#endif

#ifndef GRIDTYPE
#define GRIDTYPE Dune::YaspGrid<2>
#endif

#ifndef NAVIER_STOKES_MODEL
#define NAVIER_STOKES_MODEL NavierStokesMomentum
#endif

#ifndef DISCRETIZATION_MODEL
#define DISCRETIZATION_MODEL FaceCenteredStaggeredModel
#endif

#include <dune/grid/yaspgrid.hh>

#include <dumux/flux/fluxvariablescaching.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>

#include <dumux/discretization/fcdiamond.hh>
#include <dumux/discretization/pq1bubble.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DoneaTestMomentum { using InheritsFrom = std::tuple<NAVIER_STOKES_MODEL, DISCRETIZATION_MODEL>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DoneaTestMomentum>
{
    using type = Dumux::DoneaTestProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DoneaTestMomentum>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DoneaTestMomentum> { using type = GRIDTYPE; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DoneaTestMomentum> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::DoneaTestMomentum> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::DoneaTestMomentum> { static constexpr bool value = ENABLECACHING; };

} // end namespace Dumux::Properties

#endif
