// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ThreePThreeCTests
 * \brief The properties of the non-isothermal injection problem where water is injected into a sand
 * column with a NAPL contamination.
 */

#ifndef DUMUX_COLUMNXYLOLPROPERTIES_HH
#define DUMUX_COLUMNXYLOLPROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/h2oairxylene.hh>
#include <dumux/material/solidstates/inertsolidstate.hh>
#include <dumux/material/solidsystems/1csolid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/3p3c/model.hh>

#include "spatialparams.hh"
#include "problem.hh"

#define ISOTHERMAL 0

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct Column { using InheritsFrom = std::tuple<ThreePThreeCNI>; };
struct ColumnBox { using InheritsFrom = std::tuple<Column, BoxModel>; };
struct ColumnCCTpfa { using InheritsFrom = std::tuple<Column, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Column> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Column> { using type = ColumnProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Column>
{ using type = FluidSystems::H2OAirXylene<GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag>
struct SolidSystem<TypeTag, TTag::Column>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Component = Dumux::Components::Constant<1, Scalar>;
    using type = SolidSystems::InertSolidPhase<Scalar, Component>;
};

//! The two-phase model uses the immiscible fluid state
template<class TypeTag>
struct SolidState<TypeTag, TTag::Column>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
public:
    using type = InertSolidState<Scalar, SolidSystem>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Column>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = ColumnSpatialParams<GridGeometry, Scalar>;
};

} // end namespace Dumux::Properties

#endif
