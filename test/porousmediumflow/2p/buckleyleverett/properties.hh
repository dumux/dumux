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

namespace Dumux::Properties::TTag {

struct TwoPBuckleyLeverettTpfa
{
    using InheritsFrom = std::tuple<TwoP, CCTpfaModel>;

    using Grid = Dune::YaspGrid<2>;

    template<class TypeTag>
    using Problem = BuckleyLeverettProblem<TypeTag>;

    template<class TypeTag>
    using LocalResidual = TwoPIncompressibleLocalResidual<TypeTag>;

    using Scalar = double;

    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar>>;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<2, Scalar>>;
    using FluidSystem = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;

    template<class TypeTag>
    using SpatialParams = BuckleyLeverettSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>, Scalar>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

} // end namespace Dumux::Properties::TTag

#endif
