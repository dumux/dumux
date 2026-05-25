// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTIC_INCOMPRESSIBLE_BLOCK_PROPERTIES_HH
#define DUMUX_HYPERELASTIC_INCOMPRESSIBLE_BLOCK_PROPERTIES_HH

#include <type_traits>

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/box.hh>
#include <dumux/solidmechanics/hyperelastic/model.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct HyperelasticIncompressibleBlock
{
    using InheritsFrom = std::tuple<Hyperelastic, BoxModel>;
    using Grid = Dune::YaspGrid<3>;  // 3D cube, quarter model
    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::HyperelasticIncompressibleBlock>
{ using type = HyperelasticIncompressibleBlockProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::HyperelasticIncompressibleBlock>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = IncompressibleBlockSpatialParams<GridGeometry, Scalar>;
};

} // end namespace Dumux::Properties

#endif
