// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup StaggeredDiscretization
 *
 * \brief Defines a type tag and some properties for ree-flow models using the staggered scheme.
           This scheme features degrees of freedom at the elements' centers and intersections (faces).
 * TODO: detailed documentation and figures
 */

#ifndef DUMUX_STAGGERD_FREE_FLOW_PROPERTIES_HH
#define DUMUX_STAGGERD_FREE_FLOW_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/intersectionmapper.hh>
#include <dumux/common/defaultmappertraits.hh>

#include <dumux/discretization/staggered.hh>
#include <dumux/discretization/staggered/fvgridgeometry.hh>
#include <dumux/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/boundarytypes.hh>

#include "facevariables.hh"
#include "velocityoutput.hh"
#include "fvgridgeometrytraits.hh"
#include "gridvolumevariables.hh"

namespace Dumux {
namespace Properties {

//! Type tag for the staggered scheme specialized for free flow.
// Create new type tags
namespace TTag {
struct StaggeredFreeFlowModel { using InheritsFrom = std::tuple<StaggeredModel>; };
} // end namespace TTag

/*!
 * \brief  Set the number of equations on the faces to 1. We only consider scalar values because the velocity vector
 *         is normal to the face.
 */
template<class TypeTag>
struct NumEqFace<TypeTag, TTag::StaggeredFreeFlowModel> { static constexpr int value = 1; };

/*!
 * \brief  For free flow models, we take the number of "physical" equations
 *         (e.g. 4 for a 3D NavierStokes problem, 3 velocity components and pressure)
 *         and subtract the number of dimensions. This yields the number of equations to be
 *         solved on the cell centers. Works also for non-isothermal models.
 */
template<class TypeTag>
struct NumEqCellCenter<TypeTag, TTag::StaggeredFreeFlowModel>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    static constexpr auto dim = GridView::dimension;
    static constexpr auto numEq = ModelTraits::numEq();
public:
    static constexpr int value = numEq - dim;
};

//! The default grid geometry
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::StaggeredFreeFlowModel>
{
private:
    static constexpr auto upwindSchemeOrder = getPropValue<TypeTag, Properties::UpwindSchemeOrder>();
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Traits = StaggeredFreeFlowDefaultFVGridGeometryTraits<GridView, upwindSchemeOrder>;
public:
    using type = StaggeredFVGridGeometry<GridView, enableCache, Traits>;
};

//! The variables living on the faces
template<class TypeTag>
struct FaceVariables<TypeTag, TTag::StaggeredFreeFlowModel>
{
private:
    using FacePrimaryVariables = GetPropType<TypeTag, Properties::FacePrimaryVariables>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr auto upwindSchemeOrder = getPropValue<TypeTag, Properties::UpwindSchemeOrder>();
public:
    using type = StaggeredFaceVariables<FacePrimaryVariables, GridView::dimension, upwindSchemeOrder>;
};

//! Set the default global volume variables cache vector class
template<class TypeTag>
struct GridVolumeVariables<TypeTag, TTag::StaggeredFreeFlowModel>
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    static constexpr auto enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Traits = StaggeredGridDefaultGridVolumeVariablesTraits<Problem, VolumeVariables>;
public:
    using type = StaggeredGridVolumeVariables<Traits, enableCache>;
};

//! The velocity output
template<class TypeTag>
struct VelocityOutput<TypeTag, TTag::StaggeredFreeFlowModel>
{
    using type = StaggeredFreeFlowVelocityOutput<GetPropType<TypeTag, Properties::GridVariables>,
                                                 GetPropType<TypeTag, Properties::SolutionVector>>;
};

/*!
 * \brief  Set the order of the upwinding scheme to 1 by default.
 */
template<class TypeTag>
struct UpwindSchemeOrder<TypeTag, TTag::StaggeredFreeFlowModel> { static constexpr int value = 1; };

} // namespace Properties
} // namespace Dumux

#endif
