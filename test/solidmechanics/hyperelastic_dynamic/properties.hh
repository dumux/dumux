// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_DYNAMIC_HYPERELASTICITY_TEST_PROPERTIES_HH
#define DUMUX_DYNAMIC_HYPERELASTICITY_TEST_PROPERTIES_HH

#include <type_traits>

#include <dune/alugrid/grid.hh>
#include <dumux/discretization/box.hh>
#include <dumux/solidmechanics/hyperelastic/model.hh>

#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct DynamicHyperelasticityTest
{
    using InheritsFrom = std::tuple<Hyperelastic, BoxModel>;
    using Grid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>;
    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::DynamicHyperelasticityTest>
{ using type = DynamicHyperelasticityProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DynamicHyperelasticityTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = DefaultDynamicHyperelasticSpatialParams<GridGeometry, Scalar>;
};


} // end namespace Dumux::Properties

#endif
