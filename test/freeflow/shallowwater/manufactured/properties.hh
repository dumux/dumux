// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief The properties for the manufactured solution test
 */
#ifndef DUMUX_SHALLOWWATER_MANUFACTURED_PROPERTIES_HH
#define DUMUX_SHALLOWWATER_MANUFACTURED_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#include <dumux/discretization/cctpfa.hh>
#include <dumux/freeflow/shallowwater/model.hh>

#include "spatialparams.hh"
#include "problem.hh"

#ifndef GRIDTYPE
#define GRIDTYPE Dune::YaspGrid<2>
#endif

namespace Dumux::Properties {

namespace TTag {
struct ManufacturedSolution { using InheritsFrom = std::tuple<ShallowWater, CCTpfaModel>; };
struct ManufacturedSolutionSecondOrder { using InheritsFrom = std::tuple<ShallowWaterSecondOrder, ManufacturedSolution>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::ManufacturedSolution>
{ using type = GRIDTYPE; };

template<class TypeTag>
struct Problem<TypeTag, TTag::ManufacturedSolution>
{ using type = Dumux::ManufacturedSolutionProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::ManufacturedSolution>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
public:
    using type = ManufacturedSolutionSpatialParams<GridGeometry, Scalar, VolumeVariables>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ManufacturedSolution>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ManufacturedSolution>
{ static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
