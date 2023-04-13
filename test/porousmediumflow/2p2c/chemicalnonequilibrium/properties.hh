// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief The properties for the 2p2c chemical nonequilibrium problem.
 */

#ifndef DUMUX_TWOPTWOC_NONEQUILIBRIUM_PROPERTIES_HH
#define DUMUX_TWOPTWOC_NONEQUILIBRIUM_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/2p2c/model.hh>

#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidstates/compositional.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct TwoPTwoCChemicalNonequilibrium { using InheritsFrom = std::tuple<TwoPTwoCNINonEquil>; };
struct TwoPTwoCChemicalNonequilibriumBox { using InheritsFrom = std::tuple<TwoPTwoCChemicalNonequilibrium, BoxModel>; };
struct TwoPTwoCChemicalNonequilibriumCC { using InheritsFrom = std::tuple<TwoPTwoCChemicalNonequilibrium, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium>
{ using type = TwoPTwoCChemicalNonequilibriumProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:

    using type = FluidSystems::H2OAir<Scalar,
                                      Components::TabulatedComponent<Components::H2O<Scalar>>,
                                      FluidSystems::H2OAirDefaultPolicy</*fastButSimplifiedRelations=*/true>,
                                      true /*useKelvinEquation*/>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TwoPTwoCChemicalNonequilibriumSpatialParams<GridGeometry, Scalar>;
};

// decide which type to use for floating values (double / quad)
template<class TypeTag>
struct Scalar<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium> { using type = double; };
template<class TypeTag>
struct Formulation<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium>
{
public:
    static const TwoPFormulation value = TwoPFormulation::p0s1;
};

template<class TypeTag>
struct UseMoles<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableThermalNonEquilibrium<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium>
{ static constexpr bool value = false; };

template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium>
{ using type = FouriersLaw<TypeTag>; };

template<class TypeTag>
struct EnergyLocalResidual<TypeTag, TTag::TwoPTwoCChemicalNonequilibrium>
{ using type = Dumux::EnergyLocalResidual<TypeTag>; };

} // end namespace Dumux::Properties

#endif
