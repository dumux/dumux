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
 * \brief Defines a type tags and some fundamental grid-related properties
 */
#ifndef DUMUX_GRID_PROPERTIES_HH
#define DUMUX_GRID_PROPERTIES_HH

#include <dumux/common/properties.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dumux/common/pointsource.hh>
#include <dumux/io/gridcreator.hh>

#include <dune/common/version.hh>

namespace Dumux
{
namespace Properties
{
//! Type tag for numeric models.
NEW_TYPE_TAG(GridProperties);

//! Use the leaf grid view if not defined otherwise
SET_TYPE_PROP(GridProperties, GridView, typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView);

//! Use the DgfGridCreator by default
SET_TYPE_PROP(GridProperties, GridCreator, GridCreator<TypeTag>);

//! Use the minimal point source implementation as default
SET_TYPE_PROP(GridProperties, PointSource, PointSource<TypeTag>);

//! Use the point source helper using the bounding box tree as a default
SET_TYPE_PROP(GridProperties, PointSourceHelper, BoundingBoxTreePointSourceHelper<TypeTag>);

//! Mapper for the grid view's vertices.
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
SET_TYPE_PROP(GridProperties,
              VertexMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView)>);
#else
SET_TYPE_PROP(GridProperties,
              VertexMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGVertexLayout>);
#endif

//! Mapper for the grid view's elements.
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
SET_TYPE_PROP(GridProperties,
              ElementMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView)>);
#else
SET_TYPE_PROP(GridProperties,
              ElementMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGElementLayout>);
#endif



} // namespace Properties
} // namespace Dumux

#endif
