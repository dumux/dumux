// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup VETest
 * \brief The properties for the immiscible 2p VE test.
 */

#ifndef DUMUX_TEST_TWOPVE_PROPERTIES_HH
#define DUMUX_TEST_TWOPVE_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/ch4.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/porousmediumflow/2pve/model.hh>

#include "problem.hh"
#include "spatialparams.hh"

#ifndef GRIDTYPE // default to yasp grid if not provided by CMake
#define GRIDTYPE Dune::YaspGrid<2,Dune::EquidistantOffsetCoordinates<double,2>>
#endif

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct TwoPVEImmiscible { using InheritsFrom = std::tuple<TwoPVE>; };
struct TwoPVEImmiscibleTpfa { using InheritsFrom = std::tuple<TwoPVEImmiscible, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPVEImmiscible> { using type = GRIDTYPE; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPVEImmiscible> { using type = TwoPVETestProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPVEImmiscible>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::H2O<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePGas<Scalar, Components::CH4<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPVEImmiscible>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = TwoPTestSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct Formulation<TypeTag, TTag::TwoPVEImmiscible>
{
    static constexpr auto value = TwoPFormulation::p0s1;
};

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TwoPVEImmiscible> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TwoPVEImmiscible> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TwoPVEImmiscible> { static constexpr bool value = false; };
} // end namespace Dumux::Properties

#endif
