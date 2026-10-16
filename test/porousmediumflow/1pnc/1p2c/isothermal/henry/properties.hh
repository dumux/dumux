// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Properties of the Henry problem of Fahs et al. \cite Fahs2016
 */
#ifndef DUMUX_TEST_HENRY_FAHS_PROPERTIES_HH
#define DUMUX_TEST_HENRY_FAHS_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/flux/cvfe/boussinesqdarcyslaw.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivityconstanttortuosity.hh>
#include <dumux/material/fluidsystems/liquidphase2c.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

namespace TTag {
struct HenryFahs { using InheritsFrom = std::tuple<OnePNC, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::HenryFahs> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::HenryFahs> { using type = HenryFahsProblem<TypeTag>; };

// equal molar masses give a constant density
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::HenryFahs>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::LiquidPhaseTwoC<Scalar, Components::Constant<1, Scalar>, Components::Constant<2, Scalar>>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::HenryFahs>
{
    using type = HenryFahsSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                        GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct UseMoles<TypeTag, TTag::HenryFahs> { static constexpr bool value = false; };

// porosity times molecular diffusion coefficient (tortuosity 1)
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::HenryFahs>
{ using type = DiffusivityConstantTortuosity<GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag>
struct EnableCompositionalDispersion<TypeTag, TTag::HenryFahs> { static constexpr bool value = true; };

template<class TypeTag>
struct AdvectionType<TypeTag, TTag::HenryFahs>
{
    using type = CVFEBoussinesqDarcysLaw<GetPropType<TypeTag, Properties::Scalar>,
                                         GetPropType<TypeTag, Properties::GridGeometry>>;
};

} // end namespace Dumux::Properties

#endif
