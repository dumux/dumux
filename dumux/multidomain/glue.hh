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
 * \ingroup MultiDomain
 * \brief A class glueing two grids of potentially different dimension geometrically.
 *        Intersections are computed using axis-aligned bounding box trees
 */

#ifndef DUMUX_MULTIDOMAIN_GLUE_HH
#define DUMUX_MULTIDOMAIN_GLUE_HH

#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/intersectionentityset.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief A convenience alias for the IntersectionEntitySet of two GridViewGeometricEntitySets
 */
template<class DomainGridView, class TargetGridView, class DomainMapper, class TargetMapper>
using MultiDomainGlue = IntersectionEntitySet<GridViewGeometricEntitySet<DomainGridView, 0, DomainMapper>,
                                              GridViewGeometricEntitySet<TargetGridView, 0, TargetMapper>>;

/*!
 * \ingroup MultiDomain
 * \brief Creates the glue object containing the intersections
 *        between two grids obtained from given grid geometries.
 * \param domainGridGeometry The grid geometry of the domain
 * \param targetGridGeometry The grid geometry of the target domain
 * \return The glue object containing the intersections
 */
template<class DomainGG, class TargetGG>
MultiDomainGlue< typename DomainGG::GridView, typename TargetGG::GridView,
                 typename DomainGG::ElementMapper, typename TargetGG::ElementMapper >
makeGlue(const DomainGG& domainGridGeometry, const TargetGG& targetGridGeometry)
{
    MultiDomainGlue< typename DomainGG::GridView, typename TargetGG::GridView,
                     typename DomainGG::ElementMapper, typename TargetGG::ElementMapper > glue;
    glue.build(domainGridGeometry.boundingBoxTree(), targetGridGeometry.boundingBoxTree());
    return glue;
}

} // end namespace Dumux

#endif
