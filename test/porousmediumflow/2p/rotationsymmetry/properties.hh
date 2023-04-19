// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief The properties of 2p rotational symmetry test
 */

#ifndef DUMUX_TEST_TWOP_ROTATIONALSYMMETRY_PROPERTIES_HH
#define DUMUX_TEST_TWOP_ROTATIONALSYMMETRY_PROPERTIES_HH

#include <dune/alugrid/grid.hh>

#include <dumux/common/properties.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/extrusion.hh>

#include <dumux/porousmediumflow/2p/model.hh>

#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

namespace TTag {
struct TwoPRotationalSymmetryDome
{ using InheritsFrom = std::tuple<TwoP, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPRotationalSymmetryDome>
{ using type = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPRotationalSymmetryDome>
{ using type = TwoPRotationalSymmetryProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPRotationalSymmetryDome>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::TabulatedComponent<Components::H2O<Scalar>, false>>;
    using NonwettingPhase = FluidSystems::OnePGas<Scalar, Components::TabulatedComponent<Components::Air<Scalar>, false>>;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPRotationalSymmetryDome>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TwoPRotationalSymmetrySpatialParams<GridGeometry, Scalar>;
};

// Set the rotational symmetric grid geometry
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::TwoPRotationalSymmetryDome>
{
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    struct Traits : public BoxDefaultGridGeometryTraits<GridView> { using Extrusion = RotationalExtrusion<0>; };
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, Traits>;
};

// Caching options
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TwoPRotationalSymmetryDome> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TwoPRotationalSymmetryDome> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TwoPRotationalSymmetryDome> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
