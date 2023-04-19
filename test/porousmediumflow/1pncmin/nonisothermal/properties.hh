// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCMinTests
 * \brief Definition of the properties for thermochemical heat storage using \f$ \textnormal{CaO},
 * \textnormal{Ca} \left( \textnormal{OH} \right)_2\f$.
 */
#ifndef DUMUX_THERMOCHEM_PROBLEM_PROPERTIES_HH
#define DUMUX_THERMOCHEM_PROBLEM_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1pncmin/model.hh>

#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/components/cao2h2.hh>
#include <dumux/material/solidsystems/compositionalsolidphase.hh>

#include "problem.hh"
#include "spatialparams.hh"
#include "modifiedcao.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct ThermoChem { using InheritsFrom = std::tuple<OnePNCMinNI>; };
struct ThermoChemBox { using InheritsFrom = std::tuple<ThermoChem, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ThermoChem> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ThermoChem> { using type = ThermoChemProblem<TypeTag>; };

// The fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ThermoChem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2ON2 = FluidSystems::H2ON2<Scalar>;
    static constexpr auto phaseIdx = H2ON2::gasPhaseIdx; // simulate the air phase
    using type = FluidSystems::OnePAdapter<H2ON2, phaseIdx>;
};

// The solid system
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::ThermoChem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ComponentOne = Components::ModifiedCaO<Scalar>;
    using ComponentTwo = Components::CaO2H2<Scalar>;
    using type = SolidSystems::CompositionalSolidPhase<Scalar, ComponentOne, ComponentTwo>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::ThermoChem>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = ThermoChemSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::ThermoChem> { static constexpr bool value = true; };

} // end namespace Properties

#endif
