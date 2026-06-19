// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_1P_NETWORK_PARALLEL_PROPERTIES_HH
#define DUMUX_TEST_1P_NETWORK_PARALLEL_PROPERTIES_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

namespace TTag {
struct OnePNetwork { using InheritsFrom = std::tuple<OneP>; };
struct OnePNetworkCCTpfa { using InheritsFrom = std::tuple<OnePNetwork, CCTpfaModel>; };
struct OnePNetworkBox { using InheritsFrom = std::tuple<OnePNetwork, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::OnePNetwork> { using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::OnePNetwork> { using type = OnePNetworkProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePNetwork>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TubesTestSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePNetwork>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
};

} // end namespace Dumux::Properties

#endif
