// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPOneCTests
 * \brief Properties of the 2p1cni heat pipe test.
 */
#ifndef DUMUX_TEST_2P1C_HEATPIPE_PROPERTIES_HH
#define DUMUX_TEST_2P1C_HEATPIPE_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/2p1c/model.hh>

#include "fluidsystem.hh"
#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux {

/*!
 * \ingroup TwoPOneCTests
 * \brief Without heat conduction, the latent heat of the counter-current flow carries the whole heat flux.
 */
template<class Scalar>
struct ThermalConductivityZero
{
    template<class VolumeVariables>
    static Scalar effectiveThermalConductivity(const VolumeVariables& volVars)
    { return 0.0; }
};

} // end namespace Dumux

namespace Dumux::Properties {

namespace TTag {
struct TwoPOneCHeatPipe { using InheritsFrom = std::tuple<TwoPOneCNI>; };
struct TwoPOneCHeatPipeBox { using InheritsFrom = std::tuple<TwoPOneCHeatPipe, BoxModel>; };
struct TwoPOneCHeatPipeTpfa { using InheritsFrom = std::tuple<TwoPOneCHeatPipe, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPOneCHeatPipe> { using type = Dune::YaspGrid<1>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPOneCHeatPipe> { using type = TwoPOneCHeatPipeProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPOneCHeatPipe>
{ using type = FluidSystems::HeatPipe<GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPOneCHeatPipe>
{
    using type = TwoPOneCHeatPipeSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                               GetPropType<TypeTag, Properties::Scalar>>;
};

template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::TwoPOneCHeatPipe>
{ using type = ThermalConductivityZero<GetPropType<TypeTag, Properties::Scalar>>; };

} // end namespace Dumux::Properties

#endif
