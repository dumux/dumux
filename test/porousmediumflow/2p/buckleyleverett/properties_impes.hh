// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTests
 * \brief Properties for the pressure step of the IMPES Buckley-Leverett test.
 */
#ifndef DUMUX_TEST_TWOP_BUCKLEYLEVERETT_PROPERTIES_IMPES_HH
#define DUMUX_TEST_TWOP_BUCKLEYLEVERETT_PROPERTIES_IMPES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include "problem_impes.hh"
#include "spatialparams.hh"
#include "spatialparams_impes.hh"

namespace Dumux::Properties::TTag {

struct TwoPBuckleyLeverettImpesPressureTpfa
{
    using InheritsFrom = std::tuple<OneP, CCTpfaModel>;

    using Grid = Dune::YaspGrid<2>;

    template<class TypeTag>
    using Problem = BuckleyLeverettImpesPressureProblem<TypeTag>;

    template<class TypeTag>
    using LocalResidual = OnePIncompressibleLocalResidual<TypeTag>;

    using Scalar = double;

    // only provides the reference density and viscosity, the phase mobilities enter via the spatial params
    using FluidSystem = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar>>;

    // effective permeability K*lambda_t*mu_ref, based on the two-phase spatial params
    template<class TypeTag>
    using SpatialParams = BuckleyLeverettImpesPressureSpatialParams<
        GetPropType<TypeTag, Properties::GridGeometry>, Scalar,
        BuckleyLeverettSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>, Scalar>
    >;

    // the transmissibilities are cached and have to be recomputed after each saturation update
    // (done by calling GridVariables::init, which forces an update of the flux variables cache)
    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

} // end namespace Dumux::Properties::TTag

#endif
