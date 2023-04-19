// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \ingroup TwoPTests
 * \brief The properties for the 2p cornerpoint test.
 */

#ifndef DUMUX_TWOP_CORNERPOINT_TEST_PROPERTIES_HH
#define DUMUX_TWOP_CORNERPOINT_TEST_PROPERTIES_HH

#if HAVE_OPM_GRID
#include <opm/grid/CpGrid.hpp>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/2p/incompressiblelocalresidual.hh>

#include <dumux/material/components/trichloroethene.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct TwoPCornerPoint { using InheritsFrom = std::tuple<CCTpfaModel, TwoP>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPCornerPoint> { using type = Dune::CpGrid; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPCornerPoint> { using type = TwoPCornerPointTestProblem<TypeTag>; };

// the local residual containing the analytic derivative methods
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::TwoPCornerPoint> { using type = TwoPIncompressibleLocalResidual<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPCornerPoint>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Trichloroethene<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPCornerPoint>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = TwoPCornerPointTestSpatialParams<GridGeometry, Scalar>;
};

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TwoPCornerPoint> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TwoPCornerPoint> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TwoPCornerPoint> { static constexpr bool value = false; };
} // end namespace Properties

#else
#warning "The opm-grid module is needed to use this class!"
#endif

#endif
