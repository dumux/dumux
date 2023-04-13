// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPNCTests
 * \brief The properties of a diffusion
 *        problem which uses the isothermal 2p2c box model.
 */
#ifndef DUMUX_TWOPNC_DIFFUSION_PROPERTIES_HH
#define DUMUX_TWOPNC_DIFFUSION_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/porousmediumflow/2pnc/model.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/flux/maxwellstefanslaw.hh>

#ifndef DIFFUSIONTYPE // default to Fick's law if not set through CMake
#define DIFFUSIONTYPE FicksLaw<TypeTag>
#endif

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct TwoPNCDiffusion { using InheritsFrom = std::tuple<TwoPNC>; };
struct TwoPNCDiffusionCC { using InheritsFrom = std::tuple<TwoPNCDiffusion, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPNCDiffusion> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPNCDiffusion> { using type = TwoPNCDiffusionProblem<TypeTag>; };

// // Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPNCDiffusion>
{
    using type = FluidSystems::H2ON2<GetPropType<TypeTag, Properties::Scalar>,
                                     FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPNCDiffusion>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TwoPNCDiffusionSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::TwoPNCDiffusion> { static constexpr bool value = true; };

//! Here we set FicksLaw or TwoPNCDiffusionsLaw
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::TwoPNCDiffusion> { using type = DIFFUSIONTYPE; };

//! Set the default formulation to pw-Sn: This can be over written in the problem.
template<class TypeTag>
struct Formulation<TypeTag, TTag::TwoPNCDiffusion>
{ static constexpr auto value = TwoPFormulation::p0s1; };

} // end namespace Dumux::Properties

#endif
