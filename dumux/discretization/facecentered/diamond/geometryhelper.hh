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
template <class GridView, class ScvType, class ScvfType>
class HyperPyramidGeometryHelper
{
    using Scalar = typename GridView::ctype;
    using GlobalPosition = typename Dune::FieldVector<Scalar, GridView::dimensionworld>;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    static constexpr auto dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;

    using ScvCornerStorage = typename ScvType::Traits::CornerStorage;
    using ScvfCornerStorage = typename ScvfType::Traits::CornerStorage;

    using LocalScvIndexType = typename ScvType::Traits::LocalIndexType;
    using LocalScvfIndexType = typename ScvfType::Traits::LocalIndexType;
    using ScvPairStorage = std::vector<std::pair<LocalScvIndexType, LocalScvIndexType>>;

public:
    HyperPyramidGeometryHelper(const Element& element)
    {
        const auto elementGeom = element.geometry();
        center_ = elementGeom.center();
        const auto refElement = referenceElement(elementGeom);

        // loop over all faces
        pScv_.resize(element.subEntities(1));
        scvDofPos_.resize(element.subEntities(1));
        for(int idx=0; idx<element.subEntities(1); ++idx)
        {
            const auto numVertices = refElement.size(idx, 1, 2);
            pScv_[idx].reserve(numVertices);
            for (int localVIdx = 0; localVIdx < numVertices; ++localVIdx)
            {
                const auto vIdx = refElement.subEntity(idx, 1, localVIdx, 2);
                pScv_[idx].push_back(elementGeom.corner(vIdx));
            }
            pScv_[idx].push_back(center_);

            scvDofPos_[idx] = elementGeom.global(refElement.position(idx, 1));
        }

        // loop over all faces
        numInteriorScvf_ = element.subEntities(2);
        pScvf_.resize(numInteriorScvf_);
        for(int idx=0; idx<numInteriorScvf_; ++idx)
        {
            const auto numVertices = refElement.size(idx, 2, dim);
            pScvf_[idx].reserve(numVertices);
            for (int localVIdx = 0; localVIdx < numVertices; ++localVIdx)
            {
                const auto vIdx = refElement.subEntity(idx, 2, localVIdx, dim);
                pScvf_[idx].push_back(elementGeom.corner(vIdx));
            }
            pScvf_[idx].push_back(center_);
        }

        scvPairs_ = std::move(scvPairs(elementGeom.corners()));
    }

    //! Create a vector with the scv corners
    ScvCornerStorage getScvCorners(LocalScvIndexType localScvIdx) const
    {
        return pScv_[localScvIdx];
    }

    //! Create a vector with the scvf corners
    ScvfCornerStorage getScvfCorners(LocalScvfIndexType localScvfIdx) const
    {
        return pScvf_[localScvfIdx];
    }

    //! Create the sub control volume face geometries on the boundary
    ScvfCornerStorage getBoundaryScvfCorners(const Intersection& is) const
    {
        ScvfCornerStorage corners;
        std::copy(pScv_[is.indexInInside()].begin(), pScv_[is.indexInInside()].end()-1,
                  std::back_inserter(corners));

        return corners;
    }

    typename ScvPairStorage::value_type getFaceScvPair(LocalScvfIndexType localScvfIdx)
    {
        return scvPairs_[localScvfIdx];
    }

    LocalScvIndexType getLocalOutsideScvIdx(LocalScvIndexType localScvIdx, LocalScvfIndexType localScvfIdx)
    {
        auto scvPair = scvPairs_[localScvfIdx];
        return (scvPair.first == localScvIdx) ? scvPair.second : scvPair.first;
    }

    LocalScvfIndexType numInteriorScvf()
    {
        return numInteriorScvf_;
    }

    template<int d = dimWorld, std::enable_if_t<(d==3), int> = 0>
    GlobalPosition normal(typename ScvPairStorage::value_type scvPair, LocalScvfIndexType localScvfIdx)
    {
        const auto& p = pScvf_[localScvfIdx];
        auto normal = Dumux::crossProduct(p[1]-p[0], p[2]-p[0]);
        normal /= normal.two_norm();

        GlobalPosition v = scvDofPos_[scvPair.second]  - scvDofPos_[scvPair.first];

        const auto s = v*normal;
        if (std::signbit(s))
            normal *= -1;

        return normal;
    }

    template<int d = dimWorld, std::enable_if_t<(d==2), int> = 0>
    GlobalPosition normal(typename ScvPairStorage::value_type scvPair, LocalScvfIndexType localScvfIdx)
    {
        //! obtain normal vector by 90° counter-clockwise rotation of t
        const auto& p = pScvf_[localScvfIdx];
        const auto t = p[1] - p[0];
        GlobalPosition normal({-t[1], t[0]});
        normal /= normal.two_norm();

        GlobalPosition v = scvDofPos_[scvPair.second]  - scvDofPos_[scvPair.first];

        const auto s = v*normal;
        if (std::signbit(s))
            normal *= -1;

        return normal;
    }

private:
    //! Get pairs
    template<int d = dimWorld, std::enable_if_t<(d==3), int> = 0>
    static ScvPairStorage scvPairs(int corners)
    {
        // proceed according to number of corners of the element
        switch (corners)
        {
        case 4: // tetrahedron
        {
            return ScvPairStorage{ {0,1},
                                   {0,2},
                                   {0,3},
                                   {1,2},
                                   {1,3},
                                   {2,3} };
        }
        case 8: // hexahedron
        {
            return ScvPairStorage{ {0,2},
                                   {1,2},
                                   {0,3},
                                   {1,3},
                                   {0,4},
                                   {1,4},
                                   {2,4},
                                   {3,4},
                                   {0,5},
                                   {1,5},
                                   {2,5},
                                   {3,5} };
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Dimond scheme scv pairs for dim=" << dim
                                                                 << " dimWorld=" << dimWorld
                                                                 << " corners=" << corners);
        }
    }

    //! Get pairs
    template<int d = dimWorld, std::enable_if_t<(d==2), int> = 0>
    static ScvPairStorage scvPairs(int corners)
    {
        // proceed according to number of corners of the element
        switch (corners)
        {
        case 3: // triangle
        {
            return ScvPairStorage{ {0,1},
                                   {0,2},
                                   {1,2} };
        }
        case 4: // quadrilateral
        {
            return ScvPairStorage{ {0,2},
                                   {1,2},
                                   {0,3},
                                   {1,3} };
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Dimond scheme scv pairs for dim=" << dim
                                                                 << " dimWorld=" << dimWorld
                                                                 << " corners=" << corners);
        }
    }

    GlobalPosition center_;
    std::vector<ScvCornerStorage> pScv_;
    std::vector<ScvfCornerStorage> pScvf_;
    std::vector<GlobalPosition> scvDofPos_;
    ScvPairStorage scvPairs_;
    LocalScvfIndexType numInteriorScvf_;
};

} // end namespace Dumux

#endif
