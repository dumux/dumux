// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsTests
 * \brief The properties of a test problem for the poro-elastic model.
 */
#ifndef DUMUX_POROELASTIC_PROPERTIES_HH
#define DUMUX_POROELASTIC_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/geomechanics/poroelastic/model.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tag
namespace TTag {
struct TestPoroElastic { using InheritsFrom = std::tuple<PoroElastic, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TestPoroElastic> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TestPoroElastic> { using type = Dumux::PoroElasticProblem<TypeTag>; };

// The fluid phase consists of one constant component
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TestPoroElastic>
{
    using type = Dumux::FluidSystems::OnePLiquid< GetPropType<TypeTag, Properties::Scalar>,
                                                  Dumux::Components::Constant<0, GetPropType<TypeTag, Properties::Scalar>> >;
};

// The spatial parameters property
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TestPoroElastic>
{
    using FS = GetPropType<TypeTag, Properties::FluidSystem>;
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using type = PoroElasticSpatialParams< GetPropType<TypeTag, Properties::Scalar>,
                                           GetPropType<TypeTag, Properties::GridGeometry>,
                                           FS, PV, Indices>;
};

} // end namespace Dumux::Properties

#endif
