// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup OnePTests
 * \brief Test for the OnePModel in combination with the NI model for a conduction problem.
 *
 * The simulation domain is a tube with an elevated temperature on the left hand side.
 */

#ifndef DUMUX_1PNI_CONDUCTION_PROBLEM_PROPERTIES_HH
#define DUMUX_1PNI_CONDUCTION_PROBLEM_PROPERTIES_HH

#include <cmath>
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/1p/model.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>

#include "problem_conduction.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {
// Create new type tags
namespace TTag {
struct OnePNIConduction { using InheritsFrom = std::tuple<OnePNI>; };
struct OnePNIConductionBox { using InheritsFrom = std::tuple<OnePNIConduction, BoxModel>; };
struct OnePNIConductionCCTpfa { using InheritsFrom = std::tuple<OnePNIConduction, CCTpfaModel>; };
struct OnePNIConductionCCMpfa { using InheritsFrom = std::tuple<OnePNIConduction, CCMpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePNIConduction> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePNIConduction> { using type = OnePNIConductionProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePNIConduction>
{
    using type = FluidSystems::OnePLiquid<GetPropType<TypeTag, Properties::Scalar>,
                                          Components::H2O<GetPropType<TypeTag, Properties::Scalar>> >;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePNIConduction>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePNISpatialParams<GridGeometry, Scalar>;
};
} // end namespace Dumux

#endif
