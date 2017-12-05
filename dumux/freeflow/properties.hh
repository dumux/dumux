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
 * \ingroup Properties
 * \file
 *
 * \brief Defines a type tag and some properties for free flow models.
 */

#ifndef DUMUX_FREE_FLOW_PROPERTIES_HH
#define DUMUX_FREE_FLOW_PROPERTIES_HH

#include <dumux/discretization/staggered/freeflow/staggeredgeometryhelper.hh>
#include <dumux/discretization/staggered/freeflow/subcontrolvolumeface.hh>
#include <dumux/discretization/staggered/freeflow/facevariables.hh>
#include <dumux/implicit/staggered/primaryvariables.hh>

#include "./staggered/boundarytypes.hh"

namespace Dumux
{
namespace Properties
{
//! Type tag for models involving flow in porous media
NEW_TYPE_TAG(FreeFlow);

SET_PROP(FreeFlow, NumEq)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    static constexpr auto dim = GridView::dimension;
    static constexpr auto numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
public:
    static constexpr int value = dim + numComponents;
};

SET_INT_PROP(FreeFlow, NumEqFace, 1); //!< set the number of equations to 1
SET_INT_PROP(FreeFlow, NumEqCellCenter, 1); //!< set the number of equations to 1


//! The sub-controlvolume face
SET_PROP(FreeFlow, SubControlVolumeFace)
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
    using type = Dumux::StaggeredSubControlVolumeFace<ScvfGeometryTraits>;
};

//! The geometry helper required for the stencils, etc.
SET_PROP(FreeFlow, StaggeredGeometryHelper)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
public:
    using type = StaggeredGeometryHelper<GridView>;
};

//! The variables living on the faces
SET_TYPE_PROP(FreeFlow, FaceVariables, StaggeredFaceVariables<TypeTag>);

//! A container class used to specify values for boundary/initial conditions
SET_PROP(FreeFlow, PrimaryVariables)
{
private:
    using CellCenterBoundaryValues = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FaceBoundaryValues = Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                 GridView::dimension>;
public:
    using type = StaggeredPrimaryVariables<TypeTag, CellCenterBoundaryValues, FaceBoundaryValues>;
};

//! A container class used to specify values for sources and Neumann BCs
SET_PROP(FreeFlow, NumEqVector)
{
private:
    using CellCenterBoundaryValues = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FaceBoundaryValues = Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                 GridView::dimension>;
public:
    using type = StaggeredPrimaryVariables<TypeTag, CellCenterBoundaryValues, FaceBoundaryValues>;
};

//! Boundary types at a single degree of freedom
SET_PROP(FreeFlow, BoundaryTypes)
{
private:
    static constexpr auto size = GET_PROP_VALUE(TypeTag, NumEqCellCenter) + GET_PROP_VALUE(TypeTag, NumEqFace);
public:
    using type = StaggeredFreeFlowBoundaryTypes<size>;
};


} // namespace Properties
} // namespace Dumux

 #endif
