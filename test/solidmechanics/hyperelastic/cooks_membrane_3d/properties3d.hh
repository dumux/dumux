// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTIC_COOKS_MEMBRANE_3D_PROPERTIES_HH
#define DUMUX_HYPERELASTIC_COOKS_MEMBRANE_3D_PROPERTIES_HH

#include <dune/alugrid/grid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/solidmechanics/hyperelastic/model.hh>

#include "spatialparams.hh"
#include "problem3d.hh"

namespace Dumux::Properties {

namespace TTag {
struct CooksMembrane3D
{
    using InheritsFrom = std::tuple<Hyperelastic, BoxModel>;
    using Grid = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;
    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::CooksMembrane3D>
{ using type = CooksMembrane3DProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::CooksMembrane3D>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = CooksMembraneHyperelasticSpatialParams<GridGeometry, Scalar>;
};

} // end namespace Dumux::Properties
#endif
