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
 * \ingroup DiamondDiscretization
 * \copydoc Dumux::FaceCenteredDiamondFVGridGeometry
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_GEOMETRY_HELPER_HH

#include <array>

#include <dumux/common/math.hh>

namespace Dumux {

//! A class to create sub control volume and sub control volume face geometries per element
template <class GridView, class ScvCornerStorage, class ScvfCornerStorage>
class PyramidGeometryHelper
{
    using Scalar = typename GridView::ctype;
    using GlobalPosition = typename Dune::FieldVector<Scalar, GridView::dimensionworld>;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    static constexpr auto dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;

public:
    PyramidGeometryHelper(const Intersection& is,
                          const typename Intersection::Geometry& geometry,
                          const GlobalPosition elementCenter)
    {
        p_.reserve(geometry.corners()+1);
        for(int i=0; i<geometry.corners(); ++i)
            p_.push_back(geometry.corner(i));

        numInteriorScvf_ = p_.size()-1;

        p_.push_back(elementCenter);
        height_ = std::abs((geometry.center()-elementCenter)*is.centerUnitOuterNormal());
        area_ = geometry.volume();
    }

    //! Create a vector with the scv corners
    ScvCornerStorage getScvCorners() const
    {
        return p_;
    }

    //! Create a vector with the scvf corners
    ScvfCornerStorage getScvfCorners(unsigned int localScvfIdx) const
    {
        if constexpr (dim == 2)
            return ScvfCornerStorage{ {p_[localScvfIdx], p_.back()} };
        else if constexpr (dim == 3)
            return ScvfCornerStorage{ {p_[localScvfIdx],
                                       p_[(localScvfIdx == numInteriorScvf_-1) ?  0 : localScvfIdx+1],
                                       p_.back()} };
        else
            return ScvfCornerStorage{ {p_.back()} };
    }

    //! Create the sub control volume face geometries on the boundary
    ScvfCornerStorage getBoundaryScvfCorners() const
    {
        ScvfCornerStorage corners;
        std::copy(p_.begin(), p_.end()-1,
                  std::back_inserter(corners));

        return corners;
    }

    unsigned int numInteriorScvf()
    {
        return numInteriorScvf_;
    }

private:
    ScvCornerStorage p_; // the points needed for construction of the scv/scvf geometries
    unsigned int numInteriorScvf_;
    Scalar height_;
    Scalar area_;
};

} // end namespace Dumux

#endif
