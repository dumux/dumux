// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Property definitions for the vector-potential Boussinesq model, used
 *        by the linear-solver comparison benchmark (small grid, 2D only).
 *
 * Identical TypeTag/physics to white-noise-perturbations/properties_vorticity.hh
 * (2D branch only) -- duplicated here rather than included so this directory
 * can keep its own small-grid params.input independent of the production tests.
 */
#ifndef DUMUX_BOUSSINESQ_LINSOLVECOMPARE_VORTICITY_PROPERTIES_HH
#define DUMUX_BOUSSINESQ_LINSOLVECOMPARE_VORTICITY_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include "../vorticity/model/model.hh"
#include "../vorticity/model/modeltraits.hh"
#include "../common/problem.hh"
#include "../common/spatialparams.hh"
#include "../common/1p_boussinesq_fluidsystem.hh"

namespace Dumux::Properties {

namespace TTag {
struct BoussinesqOneSidedRB
{ using InheritsFrom = std::tuple<BoussinesqVorticityModelBox>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::BoussinesqOneSidedRB>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<double, 2>>; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::BoussinesqOneSidedRB>
{ using type = BoussinesqVorticityModelTraits<1, 2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::BoussinesqOneSidedRB>
{ using type = BoussinesqOneSidedRBProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BoussinesqOneSidedRB>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar       = GetPropType<TypeTag, Properties::Scalar>;
    using type = BoussinesqSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BoussinesqOneSidedRB>
{ using type = FluidSystems::BoussinesqFluid<GetPropType<TypeTag, Properties::Scalar>>; };

} // end namespace Dumux::Properties

#endif
