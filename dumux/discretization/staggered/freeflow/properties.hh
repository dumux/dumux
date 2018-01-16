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
#include <dumux/discretization/staggered/properties.hh>
#include <dumux/freeflow/properties.hh>

#include "subcontrolvolumeface.hh"
#include "connectivitymap.hh"
#include "facevariables.hh"
#include "boundarytypes.hh"
#include "velocityoutput.hh"

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
    static constexpr auto dim = GridView::dimension;
    static constexpr auto numEq = GET_PROP_VALUE(TypeTag, NumEq);
public:
    static constexpr int value = GET_PROP_VALUE(TypeTag, NumEq) - dim;
};

//! The default sub-controlvolume face
SET_PROP(StaggeredFreeFlowModel, SubControlVolumeFace)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;

    struct ScvfGeometryTraits
    {
        using GridIndexType = typename Grid::LeafGridView::IndexSet::IndexType;
        using LocalIndexType = unsigned int;
        using Scalar = typename Grid::ctype;
        using Geometry = typename Grid::template Codim<1>::Geometry;
        using GlobalPosition = Dune::FieldVector<Scalar, dim>;
    };

public:
    using type = FreeFlowStaggeredSubControlVolumeFace<ScvfGeometryTraits>;
};

//! The default geometry helper required for the stencils, etc.
SET_PROP(StaggeredFreeFlowModel, StaggeredGeometryHelper)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
public:
    using type = FreeFlowStaggeredGeometryHelper<GridView>;
};

//! The variables living on the faces
SET_TYPE_PROP(StaggeredFreeFlowModel, FaceVariables, StaggeredFaceVariables<TypeTag>);

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

SET_TYPE_PROP(StaggeredFreeFlowModel, AssemblyMap, StaggeredFreeFlowConnectivityMap<typename GET_PROP_TYPE(TypeTag, FVGridGeometry),
                                                                                    typename GET_PROP(TypeTag, DofTypeIndices)>);
} // namespace Properties
} // namespace Dumux

#endif
