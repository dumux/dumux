// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \brief A class glueing two distributed grids of potentially different dimension geometrically.
 *        Intersections are computed using distributed (MPI-parallel) axis-aligned bounding box trees.
 */
#ifndef DUMUX_MULTIDOMAIN_DISTRIBUTED_GLUE_HH
#define DUMUX_MULTIDOMAIN_DISTRIBUTED_GLUE_HH

#include <memory>

#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/distributedintersectionentityset.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief A convenience alias for the DistributedIntersectionEntitySet of two GridViewGeometricEntitySets
 */
template<class DomainGridView, class TargetGridView, class DomainMapper, class TargetMapper>
using DistributedMultiDomainGlue = DistributedIntersectionEntitySet<
    GridViewGeometricEntitySet<DomainGridView, 0, DomainMapper>,
    GridViewGeometricEntitySet<TargetGridView, 0, TargetMapper>
>;

/*!
 * \ingroup MultiDomain
 * \brief Creates a distributed glue object containing the intersections
 *        between two distributed grids obtained from given grid geometries.
 * \param domainGridGeometry The grid geometry of the domain
 * \param targetGridGeometry The grid geometry of the target domain
 * \return The glue object containing the intersections (distributed by domain ownership)
 * \note This is a collective operation over the processes of the grids.
 */
template<class DomainGG, class TargetGG>
DistributedMultiDomainGlue< typename DomainGG::GridView, typename TargetGG::GridView,
                            typename DomainGG::ElementMapper, typename TargetGG::ElementMapper >
makeDistributedGlue(const DomainGG& domainGridGeometry, const TargetGG& targetGridGeometry)
{
    using DomainSet = GridViewGeometricEntitySet<typename DomainGG::GridView, 0, typename DomainGG::ElementMapper>;
    using TargetSet = GridViewGeometricEntitySet<typename TargetGG::GridView, 0, typename TargetGG::ElementMapper>;

    DistributedMultiDomainGlue< typename DomainGG::GridView, typename TargetGG::GridView,
                                typename DomainGG::ElementMapper, typename TargetGG::ElementMapper > glue;
    glue.build(
        std::make_shared<DomainSet>(domainGridGeometry.gridView(), domainGridGeometry.elementMapper()),
        std::make_shared<TargetSet>(targetGridGeometry.gridView(), targetGridGeometry.elementMapper())
    );
    return glue;
}

} // end namespace Dumux

#endif
