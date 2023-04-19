// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup OnePTests
 * \brief Test for the OnePModel in combination with the NI model for a convection problem.
 *
 * The simulation domain is a tube where water with an elevated temperature is injected
 * at a constant rate on the left hand side.
 */

#ifndef DUMUX_1PNI_CONVECTION_PROBLEM_PROPERTIES_HH
#define DUMUX_1PNI_CONVECTION_PROBLEM_PROPERTIES_HH

#include <cmath>
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>

#include <dumux/porousmediumflow/1p/model.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>

#include "problem_convection.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct OnePNIConvection { using InheritsFrom = std::tuple<OnePNI>; };
struct OnePNIConvectionBox { using InheritsFrom = std::tuple<OnePNIConvection, BoxModel>; };
struct OnePNIConvectionCCTpfa { using InheritsFrom = std::tuple<OnePNIConvection, CCTpfaModel>; };
struct OnePNIConvectionCCMpfa { using InheritsFrom = std::tuple<OnePNIConvection, CCMpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePNIConvection> { using type = Dune::YaspGrid<1>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePNIConvection> { using type = OnePNIConvectionProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePNIConvection>
{
    using type = FluidSystems::OnePLiquid<GetPropType<TypeTag, Properties::Scalar>,
                                          Components::H2O<GetPropType<TypeTag, Properties::Scalar>> >;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePNIConvection>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePNISpatialParams<GridGeometry, Scalar>;
};
} // end namespace Dumux::Properties

#endif
