// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVETests
 * \brief The properties for the radial injection into a homogeneous confined aquifer.
 */

#ifndef DUMUX_TEST_TWOPVE_RADIAL_INJECTION_PROPERTIES_HH
#define DUMUX_TEST_TWOPVE_RADIAL_INJECTION_PROPERTIES_HH

#include <tuple>

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/porousmediumflow/2pve/model.hh>

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

namespace TTag {
struct TwoPVERadialInjection { using InheritsFrom = std::tuple<TwoPVE, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPVERadialInjection>
{ using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>; };

// the grid spans the radial and the vertical coordinate and is rotated about the vertical axis
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::TwoPVERadialInjection>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    struct Traits : public CCTpfaDefaultGridGeometryTraits<GridView> { using Extrusion = RotationalExtrusion<0>; };
public:
    using type = CCTpfaFVGridGeometry<GridView, enableCache, Traits>;
};

template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPVERadialInjection> { using type = TwoPVERadialInjectionProblem<TypeTag>; };

// incompressible resident and injected fluids with constant properties, the injected supercritical CO2 is represented by a gaseous phase
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPVERadialInjection>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ResidentFluid = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar>>;
    using InjectedFluid = FluidSystems::OnePGas<Scalar, Components::Constant<2, Scalar>>;
public:
    using type = FluidSystems::TwoPImmiscible<Scalar, ResidentFluid, InjectedFluid>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPVERadialInjection>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = TwoPVERadialInjectionSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct Formulation<TypeTag, TTag::TwoPVERadialInjection>
{ static constexpr auto value = TwoPFormulation::p0s1; };

// the volume variables reconstruct the mobilities of all fine-level elements in a column
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TwoPVERadialInjection> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
