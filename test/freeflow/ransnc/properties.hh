// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RANSNCTests
 * \brief The properties of the flat plate test for the multi-component staggered grid Reynolds-averaged Navier-Stokes model.
 */
#ifndef DUMUX_RANS_NC_TEST_PROPERTIES_HH
#define DUMUX_RANS_NC_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/freeflow/compositional/zeroeqncmodel.hh>
#include <dumux/freeflow/compositional/oneeqncmodel.hh>
#include <dumux/freeflow/compositional/komegancmodel.hh>
#include <dumux/freeflow/compositional/lowrekepsilonncmodel.hh>
#include <dumux/freeflow/compositional/kepsilonncmodel.hh>
#include <dumux/freeflow/compositional/sstncmodel.hh>

#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
// Base Typetag
struct RANSNCModel { using InheritsFrom = std::tuple<StaggeredFreeFlowModel>; };
// Isothermal Typetags
struct FlatPlateNCZeroEq { using InheritsFrom = std::tuple<RANSNCModel, ZeroEqNC>; };
struct FlatPlateNCOneEq { using InheritsFrom = std::tuple<RANSNCModel, OneEqNC>; };
struct FlatPlateNCKOmega { using InheritsFrom = std::tuple<RANSNCModel, KOmegaNC>; };
struct FlatPlateNCLowReKEpsilon { using InheritsFrom = std::tuple<RANSNCModel, LowReKEpsilonNC>; };
struct FlatPlateNCKEpsilon { using InheritsFrom = std::tuple<RANSNCModel, KEpsilonNC>; };
struct FlatPlateNCSST { using InheritsFrom = std::tuple<RANSNCModel, SSTNC>; };
// Isothermal Typetags
struct FlatPlateNCNIZeroEq { using InheritsFrom = std::tuple<RANSNCModel, ZeroEqNCNI>; };
struct FlatPlateNCNIOneEq { using InheritsFrom = std::tuple<RANSNCModel, OneEqNCNI>; };
struct FlatPlateNCNIKOmega { using InheritsFrom = std::tuple<RANSNCModel, KOmegaNCNI>; };
struct FlatPlateNCNILowReKEpsilon { using InheritsFrom = std::tuple<RANSNCModel, LowReKEpsilonNCNI>; };
struct FlatPlateNCNIKEpsilon { using InheritsFrom = std::tuple<RANSNCModel, KEpsilonNCNI>; };
struct FlatPlateNCNISST { using InheritsFrom = std::tuple<RANSNCModel, SSTNCNI>; };
} // end namespace TTag

// The fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RANSNCModel>
{
  using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
  static constexpr auto phaseIdx = H2OAir::gasPhaseIdx; // simulate the air phase
  using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

// replace the main component balance eq with a total balance eq
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::RANSNCModel> { static constexpr int value = 0; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::RANSNCModel>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::RANSNCModel> { using type = Dumux::FlatPlateNCTestProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::RANSNCModel> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::RANSNCModel> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::RANSNCModel> { static constexpr bool value = true; };

template<class TypeTag>
struct UseMoles<TypeTag, TTag::RANSNCModel> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
