// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Volume of a network of tubes embedded in a bulk grid, per bulk element
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_LOWDIMVOLUME_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_LOWDIMVOLUME_HH

#include <cmath>
#include <cstddef>
#include <type_traits>
#include <vector>

#include <dumux/multidomain/glue.hh>

namespace Dumux::EmbeddedCoupling {

/*!
 * \ingroup EmbeddedCoupling
 * \brief The volume that a network of tubes occupies in each bulk element
 *
 * The network is the one-dimensional low-dimensional grid. Each element of it is the centerline
 * of a straight tube with circular cross-section and the radius given by the radius function.
 * Each intersection of the glue is a piece of a centerline inside the bulk domain. The volume
 * \f$ \pi r^2 l \f$ of the cylinder around the piece of length \f$ l \f$ is split equally among
 * the bulk elements that contain the piece (more than one if the piece lies on a face, edge or
 * vertex of the bulk grid), so the sum over all bulk elements is the sum of the cylinder volumes.
 *
 * The result is an approximation that is wrong in the following cases:
 * - The cylinder is attributed only to the bulk elements that contain the centerline. The volume
 *   of the cylinder that lies in other bulk elements is attributed to the wrong elements, so the
 *   result is only accurate if the radius is small compared to the size of the bulk elements.
 * - The equal split among bulk elements that share the centerline ignores the actual shares of
 *   the cylinder volume in them.
 * - The volume of a part of a tube that lies outside of the bulk domain is not counted.
 * - The volume where the cylinders of several elements overlap, e.g. at junctions, bifurcations
 *   and bends of the network, is counted once for each of the elements.
 * - The ends of each cylinder are flat. The cylinder has the constant radius of its element,
 *   a radius that varies along the network is approximated by one value per element.
 *
 * \param glue the intersections of the bulk grid and the low-dimensional grid
 * \param bulkGridGeometry the grid geometry of the bulk domain
 * \param lowDimGridGeometry the grid geometry of the low-dimensional domain (one-dimensional grid)
 * \param radius a function that returns the tube radius for the index of a low-dimensional element
 * \return the volume for each bulk element, indexed with the element mapper of the bulk grid geometry
 */
template<class Glue, class BulkGridGeometry, class LowDimGridGeometry, class RadiusFunction>
auto lowDimVolumePerBulkElement(const Glue& glue,
                                const BulkGridGeometry& bulkGridGeometry,
                                const LowDimGridGeometry& lowDimGridGeometry,
                                const RadiusFunction& radius)
{
    static_assert(LowDimGridGeometry::GridView::dimension == 1,
                  "The low-dimensional volume is defined for networks of tubes (one-dimensional grids)");

    using Scalar = std::decay_t<decltype(radius(std::size_t{}))>;
    std::vector<Scalar> volume(bulkGridGeometry.gridView().size(0), 0.0);

    for (const auto& is : intersections(glue))
    {
        // all target entities of an intersection are identical
        const auto lowDimElementIdx = lowDimGridGeometry.elementMapper().index(is.targetEntity(0));
        const auto r = radius(lowDimElementIdx);
        const Scalar cylinderVolume = M_PI*r*r*is.geometry().volume();

        for (std::size_t i = 0; i < is.numDomainNeighbors(); ++i)
        {
            const auto bulkElementIdx = bulkGridGeometry.elementMapper().index(is.domainEntity(i));
            volume[bulkElementIdx] += cylinderVolume/is.numDomainNeighbors();
        }
    }

    return volume;
}

/*!
 * \ingroup EmbeddedCoupling
 * \brief The volume that a network of tubes occupies in each bulk element
 * \note This overload intersects the two grids. See the overload taking the glue for the
 *       description of the result and its limitations.
 */
template<class BulkGridGeometry, class LowDimGridGeometry, class RadiusFunction>
auto lowDimVolumePerBulkElement(const BulkGridGeometry& bulkGridGeometry,
                                const LowDimGridGeometry& lowDimGridGeometry,
                                const RadiusFunction& radius)
{
    return lowDimVolumePerBulkElement(makeGlue(bulkGridGeometry, lowDimGridGeometry),
                                      bulkGridGeometry, lowDimGridGeometry, radius);
}

} // end namespace Dumux::EmbeddedCoupling

#endif
