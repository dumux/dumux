// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
