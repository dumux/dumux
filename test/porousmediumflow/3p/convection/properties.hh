// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ThreePTests
 * \brief The properties of the test for the ThreePModel in combination with the NI model for a convection problem.
 */

#ifndef DUMUX_3PNI_CONVECTION_PROPERTIES_HH
#define DUMUX_3PNI_CONVECTION_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/3p/model.hh>
#include <dumux/material/fluidsystems/h2oairmesitylene.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidmatrixinteractions/3p/thermalconductivitysomerton3p.hh>

#include "../conduction/spatialparams.hh"  // reuse the conduction spatialParams
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct ThreePNIConvection { using InheritsFrom = std::tuple<ThreePNI>; };
struct ThreePNIConvectionBox { using InheritsFrom = std::tuple<ThreePNIConvection, BoxModel>; };
struct ThreePNIConvectionCCTpfa { using InheritsFrom = std::tuple<ThreePNIConvection, CCTpfaModel>; };
struct ThreePNIConvectionCCMpfa { using InheritsFrom = std::tuple<ThreePNIConvection, CCMpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ThreePNIConvection> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ThreePNIConvection> { using type = ThreePNIConvectionProblem<TypeTag>; };


// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ThreePNIConvection>
{ using type = FluidSystems::H2OAirMesitylene<GetPropType<TypeTag, Properties::Scalar>>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::ThreePNIConvection>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = ThreePNISpatialParams<GridGeometry, Scalar>;
};

} // end namespace Dumux::Properties

#endif
