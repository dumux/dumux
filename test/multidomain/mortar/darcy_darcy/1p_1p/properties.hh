// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \ingroup MultiDomainTests
 * \brief The properties for the single-phase darcy-darcy mortar-coupling test
 */
#ifndef DUMUX_MORTAR_DARCY_ONEP_TEST_PROPERTIES_HH
#define DUMUX_MORTAR_DARCY_ONEP_TEST_PROPERTIES_HH

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/1p/model.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "problem.hh"
#include "spatialparams.hh"
#include "property_declarations.hh"

namespace Dumux::Properties {

namespace TTag {
struct OnePDarcyMortar { using InheritsFrom = std::tuple<OneP>; };
struct OnePDarcyMortarTpfa { using InheritsFrom = std::tuple<OnePDarcyMortar, CCTpfaModel>; };
struct OnePDarcyMortarBox { using InheritsFrom = std::tuple<OnePDarcyMortar, BoxModel>; };
} // namespace TTag

template<class TypeTag>
struct MortarGrid<TypeTag, TTag::OnePDarcyMortar> { using type = Dune::FoamGrid<1, 2>; };

template<class TypeTag>
struct MortarSolutionVector<TypeTag, TTag::OnePDarcyMortar> { using type = Dune::BlockVector<Dune::FieldVector<double, 1>>; };

template<class TypeTag>
struct Grid<TypeTag, TTag::OnePDarcyMortar> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::OnePDarcyMortar> { using type = DarcyProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePDarcyMortar>
{
    using type = DarcySpatialParams<
        GetPropType<TypeTag, Properties::GridGeometry>,
        GetPropType<TypeTag, Properties::Scalar>
    >;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePDarcyMortar>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<0, Scalar> >;
};

}  // namespace Dumux::Properties

#endif
