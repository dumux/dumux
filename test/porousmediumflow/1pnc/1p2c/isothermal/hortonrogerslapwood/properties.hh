// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Properties of the Horton-Rogers-Lapwood benchmark
 */
#ifndef DUMUX_TEST_HORTON_ROGERS_LAPWOOD_PROPERTIES_HH
#define DUMUX_TEST_HORTON_ROGERS_LAPWOOD_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/flux/cvfe/boussinesqdarcyslaw.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh>
#include <dumux/material/fluidsystems/liquidphase2c.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/fvspatialparams1pconstant.hh>

#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct HortonRogersLapwood { using InheritsFrom = std::tuple<OnePNC, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::HortonRogersLapwood> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::HortonRogersLapwood> { using type = HortonRogersLapwoodProblem<TypeTag>; };

// components with equal molar masses give a constant (reference) density
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::HortonRogersLapwood>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::LiquidPhaseTwoC<Scalar, Components::Constant<1, Scalar>, Components::Constant<2, Scalar>>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::HortonRogersLapwood>
{
    using type = FVPorousMediumFlowSpatialParamsOnePConstant<GetPropType<TypeTag, Properties::GridGeometry>,
                                                             GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct UseMoles<TypeTag, TTag::HortonRogersLapwood> { static constexpr bool value = false; };

// the total mass balance keeps the boundaries with fixed mass fraction impermeable
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::HortonRogersLapwood> { static constexpr int value = 0; };

template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::HortonRogersLapwood>
{ using type = DiffusivityConstantTortuosity<GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag>
struct AdvectionType<TypeTag, TTag::HortonRogersLapwood>
{
    using type = CVFEBoussinesqDarcysLaw<GetPropType<TypeTag, Properties::Scalar>,
                                         GetPropType<TypeTag, Properties::GridGeometry>>;
};

} // end namespace Dumux::Properties

#endif
