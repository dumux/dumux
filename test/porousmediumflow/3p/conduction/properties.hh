// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ThreePTests
 * \brief The properties of the 3pni conduction problem.
 */

#ifndef DUMUX_3PNI_CONDUCTION_PROPERTIES_HH
#define DUMUX_3PNI_CONDUCTION_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/3p/model.hh>
#include <dumux/material/fluidsystems/h2oairmesitylene.hh>

#include <dumux/material/fluidmatrixinteractions/3p/thermalconductivitysomerton3p.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct ThreePNIConduction { using InheritsFrom = std::tuple<ThreePNI>; };
struct ThreePNIConductionBox { using InheritsFrom = std::tuple<ThreePNIConduction, BoxModel>; };
struct ThreePNIConductionCCTpfa { using InheritsFrom = std::tuple<ThreePNIConduction, CCTpfaModel>; };
struct ThreePNIConductionCCMpfa { using InheritsFrom = std::tuple<ThreePNIConduction, CCMpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ThreePNIConduction> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ThreePNIConduction> { using type = ThreePNIConductionProblem<TypeTag>; };


// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ThreePNIConduction>
{ using type = FluidSystems::H2OAirMesitylene<GetPropType<TypeTag, Properties::Scalar>>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::ThreePNIConduction>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = ThreePNISpatialParams<GridGeometry, Scalar>;
};

}// end namespace Dumux::Properties

#endif
