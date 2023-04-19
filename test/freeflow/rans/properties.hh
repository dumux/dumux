// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RANSTests
 * \brief Pipe flow test for the staggered grid RANS model
 *
 * This test simulates pipe flow experiments performed by John Laufer in 1954
 * \cite Laufer1954a.
 */
#ifndef DUMUX_PIPE_LAUFER_PROPERTIES_HH
#define DUMUX_PIPE_LAUFER_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/turbulenceproperties.hh>
#include <dumux/freeflow/turbulencemodel.hh>

#include <dumux/freeflow/rans/problem.hh>
#include <dumux/freeflow/rans/zeroeq/model.hh>
#include <dumux/freeflow/rans/oneeq/problem.hh>
#include <dumux/freeflow/rans/oneeq/model.hh>
#include <dumux/freeflow/rans/twoeq/kepsilon/model.hh>
#include <dumux/freeflow/rans/twoeq/kepsilon/problem.hh>
#include <dumux/freeflow/rans/twoeq/komega/model.hh>
#include <dumux/freeflow/rans/twoeq/komega/problem.hh>
#include <dumux/freeflow/rans/twoeq/sst/model.hh>
#include <dumux/freeflow/rans/twoeq/sst/problem.hh>
#include <dumux/freeflow/rans/twoeq/lowrekepsilon/problem.hh>
#include <dumux/freeflow/rans/twoeq/lowrekepsilon/model.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
// Base Typetag
struct RANSModel { using InheritsFrom = std::tuple<StaggeredFreeFlowModel>; };
// Isothermal Typetags
struct PipeLauferZeroEq { using InheritsFrom = std::tuple<RANSModel, ZeroEq>; };
struct PipeLauferOneEq { using InheritsFrom = std::tuple<RANSModel, OneEq>; };
struct PipeLauferKOmega { using InheritsFrom = std::tuple<RANSModel, KOmega>; };
struct PipeLauferLowReKEpsilon { using InheritsFrom = std::tuple<RANSModel, LowReKEpsilon>; };
struct PipeLauferKEpsilon { using InheritsFrom = std::tuple<RANSModel, KEpsilon>; };
struct PipeLauferSST { using InheritsFrom = std::tuple<RANSModel, SST>; };
// Non-Isothermal Typetags
struct PipeLauferNIZeroEq { using InheritsFrom = std::tuple<RANSModel, ZeroEqNI>; };
struct PipeLauferNIOneEq { using InheritsFrom = std::tuple<RANSModel, OneEqNI>; };
struct PipeLauferNIKOmega { using InheritsFrom = std::tuple<RANSModel, KOmegaNI>; };
struct PipeLauferNILowReKEpsilon { using InheritsFrom = std::tuple<RANSModel, LowReKEpsilonNI>; };
struct PipeLauferNIKEpsilon { using InheritsFrom = std::tuple<RANSModel, KEpsilonNI>; };
struct PipeLauferNISST { using InheritsFrom = std::tuple<RANSModel, SSTNI>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RANSModel>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePGas<Scalar, Components::Air<Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::RANSModel>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::RANSModel>
{ using type = Dumux::PipeLauferProblem<TypeTag>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::RANSModel> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::RANSModel> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::RANSModel> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
