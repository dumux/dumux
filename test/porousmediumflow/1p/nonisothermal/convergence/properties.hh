// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Properties of the convergence test of the non-isothermal single-phase model
 */
#ifndef DUMUX_ONEPNI_CONVERGENCE_PROPERTIES_HH
#define DUMUX_ONEPNI_CONVERGENCE_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/pq3.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/porousmediumflow/1p/model.hh>

#include "spatialparams.hh"
#include "problem.hh"

#ifndef GRIDTYPE
#define GRIDTYPE Dune::YaspGrid<2>
#endif

namespace Dumux::Properties {

namespace TTag {
struct OnePNIConvergence { using InheritsFrom = std::tuple<OnePNI>; };
struct OnePNIConvergenceBox { using InheritsFrom = std::tuple<OnePNIConvergence, BoxModel>; };
struct OnePNIConvergencePQ1Bubble { using InheritsFrom = std::tuple<OnePNIConvergence, PQ1BubbleHybridModel>; };
struct OnePNIConvergencePQ2 { using InheritsFrom = std::tuple<OnePNIConvergence, PQ2HybridModel>; };
struct OnePNIConvergencePQ3 { using InheritsFrom = std::tuple<OnePNIConvergence, PQ3HybridModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::OnePNIConvergence> { using type = GRIDTYPE; };

template<class TypeTag>
struct Problem<TypeTag, TTag::OnePNIConvergence> { using type = OnePNIConvergenceProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePNIConvergence>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePNIConvergence>
{
    using type = OnePNIConvergenceSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                                GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::OnePNIConvergence> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
