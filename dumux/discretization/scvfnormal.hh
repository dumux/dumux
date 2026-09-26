// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Normal vector of a sub-control volume face of an element-local dual grid
 */
#ifndef DUMUX_DISCRETIZATION_SCVF_NORMAL_HH
#define DUMUX_DISCRETIZATION_SCVF_NORMAL_HH

#include <dumux/common/math.hh>
#include <dumux/geometry/normal.hh>

namespace Dumux::Detail {

/*!
 * \ingroup Discretization
 * \brief Unit normal of an interior sub-control volume face, up to sign
 * \param geo the geometry of the element the sub-control volume face belongs to
 * \param corners the corners of the sub-control volume face
 *
 * The face is a codimension-one sub-entity of the element-local dual grid, so its
 * normal is normal to the face and tangential to the element. For a surface grid
 * (dim < dimWorld) the second condition is a restriction: the normal is not the
 * surface normal but lies in the plane of the element. Callers orient the result
 * themselves, since the point pair defining the sign differs between schemes.
 */
template<class ElementGeometry, class CornerStorage>
inline auto scvfUnitNormal(const ElementGeometry& geo, const CornerStorage& corners)
{
    using GlobalPosition = typename ElementGeometry::GlobalCoordinate;
    static constexpr int dim = ElementGeometry::mydimension;
    static constexpr int dimWorld = ElementGeometry::coorddimension;
    static_assert(dim > 1 && dim <= 3 && dimWorld <= 3);

    const auto normalize = [](GlobalPosition n){ n /= n.two_norm(); return n; };

    if constexpr (dim == 2 && dimWorld == 2)
        return normalize(Dumux::normal(GlobalPosition(corners[1] - corners[0])));

    else if constexpr (dim == 2 && dimWorld == 3)
    {
        // A segment's two corners do not determine the element plane.
        const auto elementNormal = crossProduct(
            GlobalPosition(geo.corner(1) - geo.corner(0)),
            GlobalPosition(geo.corner(2) - geo.corner(0))
        );
        return normalize(crossProduct(elementNormal, GlobalPosition(corners[1] - corners[0])));
    }

    else
        return normalize(crossProduct(
            GlobalPosition(corners[1] - corners[0]),
            GlobalPosition(corners[2] - corners[0])
        ));
}

} // end namespace Dumux::Detail

#endif
