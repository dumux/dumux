// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Properties of the hydrostatic one-phase test
 */

#ifndef DUMUX_TEST_ONEP_HYDROSTATIC_PROPERTIES_HH
#define DUMUX_TEST_ONEP_HYDROSTATIC_PROPERTIES_HH

#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#endif

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct OnePHydrostatic { using InheritsFrom = std::tuple<OneP>; };
struct OnePHydrostaticTpfa { using InheritsFrom = std::tuple<OnePHydrostatic, CCTpfaModel>; };
struct OnePHydrostaticBox { using InheritsFrom = std::tuple<OnePHydrostatic, BoxModel>; };
struct OnePHydrostaticSurfaceTpfa { using InheritsFrom = std::tuple<OnePHydrostatic, CCTpfaModel>; };
struct OnePHydrostaticNetworkTpfa { using InheritsFrom = std::tuple<OnePHydrostatic, CCTpfaModel>; };
} // end namespace TTag

#if HAVE_DUNE_UGGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePHydrostatic>
{ using type = Dune::UGGrid<2>; };
#endif

#if HAVE_DUNE_FOAMGRID
// two-dimensional surface grid embedded in 3D
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePHydrostaticSurfaceTpfa>
{ using type = Dune::FoamGrid<2, 3>; };

// one-dimensional network embedded in 2D
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePHydrostaticNetworkTpfa>
{ using type = Dune::FoamGrid<1, 2>; };
#endif

template<class TypeTag>
struct Problem<TypeTag, TTag::OnePHydrostatic>
{ using type = OnePHydrostaticProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePHydrostatic>
{
    using type = OnePHydrostaticSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                              GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePHydrostatic>
{ using type = OnePIncompressibleLocalResidual<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePHydrostatic>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar>>;
};

} // end namespace Dumux::Properties

#endif
