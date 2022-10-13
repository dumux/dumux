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
 * \ingroup CCTpfaDiscretization
 * \brief Helper functionality to get outer normal vectors on element facets.
 */
#ifndef DUMUX_GEOMETRY_NORMAL_HH
#define DUMUX_GEOMETRY_NORMAL_HH

#include <algorithm>

#include <dune/geometry/referenceelement.hh>

namespace Dumux {

// TODO: proper helper (with compile-time overloads)
//       this should also get rid of the currently appearing compiler warnings
template<typename Geometry>
auto getNormal(const Geometry& geo, unsigned int facetIdx)
{
    static_assert(Geometry::coorddimension == 2);
    const auto scale = [] (auto&& v) {
        v /= v.two_norm();
        return v;
    };
    const auto rotateClockWiseAndScale = [&] (auto&& v) {
        std::swap(v[0], v[1]);
        v[1] *= -1.0;
        scale(v);
        return v;
    };

    using Dune::referenceElement;
    const auto refElement = referenceElement(geo);
    if (geo.type() == Dune::GeometryTypes::quadrilateral)
    {
        switch (facetIdx)
        {
            case 0: return rotateClockWiseAndScale(
                geo.global(refElement.position(0, 2)) -
                geo.global(refElement.position(2, 2))
            );
            case 1: return rotateClockWiseAndScale(
                geo.global(refElement.position(3, 2)) -
                geo.global(refElement.position(1, 2))
            );
            case 2: return rotateClockWiseAndScale(
                geo.global(refElement.position(1, 2)) -
                geo.global(refElement.position(0, 2))
            );
            case 3: return rotateClockWiseAndScale(
                geo.global(refElement.position(2, 2)) -
                geo.global(refElement.position(3, 2))
            );
        }
    }
    if (geo.type() == Dune::GeometryTypes::line)
        return scale(geo.global(refElement.position(facetIdx, 1)) - geo.center());

    DUNE_THROW(Dune::NotImplemented, "TODO");
}

} // end namespace Dumux

#endif
