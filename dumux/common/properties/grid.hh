// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \file
 * \ingroup Properties
 * \brief Defines a type tags and some fundamental grid-related properties
 */
#ifndef DUMUX_GRID_PROPERTIES_HH
#define DUMUX_GRID_PROPERTIES_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/pointsource.hh>

namespace Dumux {
namespace Properties {

namespace TTag {
//! Type tag for numeric models.
struct GridProperties {};
}

//! Use the minimal point source implementation as default
template<class TypeTag>
struct PointSource<TypeTag, TTag::GridProperties>
{
private:
    using SourceValues = GetPropType<TypeTag, Properties::NumEqVector>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using GlobalPosition = typename Dune::FieldVector<typename GridView::ctype, GridView::dimensionworld>;
public:
    using type = Dumux::PointSource<GlobalPosition, SourceValues>;
};

//! Use the point source helper using the bounding box tree as a default
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::GridProperties> { using type = BoundingBoxTreePointSourceHelper; };

} // namespace Properties
} // namespace Dumux

#endif
