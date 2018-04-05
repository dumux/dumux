// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \ingroup StaggeredDiscretization
 * \file
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

#include <dumux/discretization/staggered/properties.hh>
#include <dumux/discretization/staggered/fvgridgeometry.hh>
#include <dumux/freeflow/properties.hh>

#include "facevariables.hh"
#include "boundarytypes.hh"
#include "velocityoutput.hh"
#include "fvgridgeometrytraits.hh"

namespace Dumux
{

namespace Properties
{

//! Type tag for the staggered scheme specialized for free flow.
NEW_TYPE_TAG(StaggeredFreeFlowModel, INHERITS_FROM(StaggeredModel));

/*!
 * \brief  Set the number of equations on the faces to 1. We only consider scalar values because the velocity vector
 *         is normal to the face.
 */
SET_INT_PROP(StaggeredFreeFlowModel, NumEqFace, 1);

/*!
 * \brief  For free flow models, we take the number of "physical" equations
 *         (e.g. 4 for a 3D NavierStokes problem, 3 velocity components and pressure)
 *         and substract the number of dimensions. This yields the number of equations to be
 *         solved on the cell centers. Works also for non-isothermal models.
 */
SET_PROP(StaggeredFreeFlowModel, NumEqCellCenter)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    static constexpr auto dim = GridView::dimension;
    static constexpr auto numEq = ModelTraits::numEq();
public:
    static constexpr int value = numEq - dim;
};

//! The default fv grid geometry
SET_PROP(StaggeredFreeFlowModel, FVGridGeometry)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Traits = StaggeredFreeFlowDefaultFVGridGeometryTraits<GridView>;
    static constexpr bool enableCache = GET_PROP_VALUE(TypeTag, EnableFVGridGeometryCache);
public:
    using type = StaggeredFVGridGeometry<GridView, enableCache, Traits>;
};

//! The variables living on the faces
SET_PROP(StaggeredFreeFlowModel, FaceVariables)
{
private:
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
public:
    using type = StaggeredFaceVariables<FacePrimaryVariables, GridView::dimension>;
};

//! Boundary types at a single degree of freedom
SET_PROP(StaggeredFreeFlowModel, BoundaryTypes)
{
private:
    static constexpr auto size = GET_PROP_VALUE(TypeTag, NumEqCellCenter) + GET_PROP_VALUE(TypeTag, NumEqFace);
public:
    using type = StaggeredFreeFlowBoundaryTypes<size>;
};

//! The velocity output
SET_TYPE_PROP(StaggeredFreeFlowModel, VelocityOutput, StaggeredFreeFlowVelocityOutput<TypeTag>);

} // namespace Properties
} // namespace Dumux

#endif
