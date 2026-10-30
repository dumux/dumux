// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief The properties of the non-isothermal single-phase test with two layers of different thermal conductivity
 */
#ifndef DUMUX_ONEPNI_LAYERED_PROPERTIES_HH
#define DUMUX_ONEPNI_LAYERED_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/porousmediumflow/1p/model.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct OnePNILayered { using InheritsFrom = std::tuple<OnePNI>; };
struct OnePNILayeredBox { using InheritsFrom = std::tuple<OnePNILayered, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::OnePNILayered> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::OnePNILayered> { using type = OnePNILayeredProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePNILayered>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePNILayered>
{
    using type = OnePNILayeredSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                            GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::OnePNILayered> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
