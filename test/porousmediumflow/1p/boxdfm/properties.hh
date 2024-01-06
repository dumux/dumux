// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \ingroup TwoPTests
 * \brief The properties for the incompressible 1p-boxdfm test.
 */

#ifndef DUMUX_INCOMPRESSIBLE_ONEPBOXDFM_TEST_PROPERTIES_HH
#define DUMUX_INCOMPRESSIBLE_ONEPBOXDFM_TEST_PROPERTIES_HH

#include <dune/alugrid/grid.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/boxdfm/model.hh>

#include "problem.hh"
#include "spatialparams.hh"

#ifndef GRIDGEOMETRYCACHING
#define GRIDGEOMETRYCACHING 0
#endif

namespace Dumux::Properties {

// we need to derive first from twop and then from the box-dfm Model
// because the flux variables cache type of TwoP is overwritten in BoxDfmModel
// Create new type tags
namespace TTag {
struct OnePIncompressibleBoxDfm { using InheritsFrom = std::tuple<BoxDfmModel, OneP>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePIncompressibleBoxDfm> { using type = Dune::ALUGrid<2, 2, Dune::simplex , Dune::conforming>; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePIncompressibleBoxDfm> { using type = OnePTestProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePIncompressibleBoxDfm>
{
private:
    using FVG = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = OnePTestSpatialParams<FVG, Scalar>;
};

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePIncompressibleBoxDfm>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar>>;
};

// Caching preferences
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePIncompressibleBoxDfm> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePIncompressibleBoxDfm> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePIncompressibleBoxDfm> { static constexpr bool value = GRIDGEOMETRYCACHING; };

} // end namespace Dumux::Properties

#endif
