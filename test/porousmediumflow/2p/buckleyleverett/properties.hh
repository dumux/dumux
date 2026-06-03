// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTests
 * \brief Properties for the Buckley-Leverett two-phase test.
 */
#ifndef DUMUX_TEST_TWOP_BUCKLEYLEVERETT_PROPERTIES_HH
#define DUMUX_TEST_TWOP_BUCKLEYLEVERETT_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/2p/incompressiblelocalresidual.hh>

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

namespace TTag {
struct TwoPBuckleyLeverett { using InheritsFrom = std::tuple<TwoP>; };
struct TwoPBuckleyLeverettTpfa { using InheritsFrom = std::tuple<TwoPBuckleyLeverett, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPBuckleyLeverett> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPBuckleyLeverett> { using type = BuckleyLeverettProblem<TypeTag>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::TwoPBuckleyLeverett> { using type = TwoPIncompressibleLocalResidual<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPBuckleyLeverett>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar>>;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<2, Scalar>>;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPBuckleyLeverett>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = BuckleyLeverettSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TwoPBuckleyLeverett> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TwoPBuckleyLeverett> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TwoPBuckleyLeverett> { static constexpr bool value = false; };

} // end namespace Dumux::Properties

#endif
