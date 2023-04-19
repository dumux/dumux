// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ThreePThreeCTests
 * \brief The properties of the isothermal NAPL infiltration problem: LNAPL contaminates
 *        the unsaturated and the saturated groundwater zone.
 */

#ifndef DUMUX_INFILTRATION_THREEPTHREEC_PROPERTIES_HH
#define DUMUX_INFILTRATION_THREEPTHREEC_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/3p3c/model.hh>
#include <dumux/material/fluidsystems/h2oairmesitylene.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct InfiltrationThreePThreeC { using InheritsFrom = std::tuple<ThreePThreeC>; };
struct InfiltrationThreePThreeCBox { using InheritsFrom = std::tuple<InfiltrationThreePThreeC, BoxModel>; };
struct InfiltrationThreePThreeCCCTpfa { using InheritsFrom = std::tuple<InfiltrationThreePThreeC, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::InfiltrationThreePThreeC> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::InfiltrationThreePThreeC> { using type = InfiltrationThreePThreeCProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::InfiltrationThreePThreeC>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = InfiltrationThreePThreeCSpatialParams<GridGeometry, Scalar>;
};

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::InfiltrationThreePThreeC>
{ using type = FluidSystems::H2OAirMesitylene<GetPropType<TypeTag, Properties::Scalar>>; };

} // end namespace Dumux::Properties

#endif
