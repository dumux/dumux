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
 * \ingroup BoxDiscretization
 * \brief Helper class constructing the dual grid finite volume geometries
 *        for the box discretizazion method
 */
#ifndef DUMUX_DISCRETIZATION_BOX_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_BOX_GEOMETRY_HELPER_HH

#include <array>

#include <dumux/common/math.hh>

namespace Dumux {

//! Create sub control volumes and sub control volume face geometries
template<class GridView, int dim, class ScvType, class ScvfType>
class BoxGeometryHelper;

//! A class to create sub control volume and sub control volume face geometries per element
template <class GridView, class ScvType, class ScvfType>
class BoxGeometryHelper<GridView, 1, ScvType, ScvfType>
{
private:
    using Scalar = typename GridView::ctype;
    using GlobalPosition = typename Dune::FieldVector<Scalar, GridView::dimensionworld>;
    using ScvCornerStorage = typename ScvType::Traits::CornerStorage;
    using ScvfCornerStorage = typename ScvfType::Traits::CornerStorage;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    //! the maximum number of helper points used to construct the geometries
    //! Using a statically sized point array is much faster than dynamic allocation
    static constexpr int maxPoints = 3;

public:

    BoxGeometryHelper(const typename Element::Geometry& geometry)
    : elementGeometry_(geometry), corners_(geometry.corners()),
      p_({geometry.center(), geometry.corner(0), geometry.corner(1)}) {}

    //! Create a vector with the scv corners
    ScvCornerStorage getScvCorners(unsigned int localScvIdx) const
    {
        //! Only build the maps the first time we call this function
        static const std::uint8_t map[2][2] =
        {
            {1, 0},
            {2, 0}
        };

        return ScvCornerStorage{ {p_[map[localScvIdx][0]],
                                  p_[map[localScvIdx][1]]} };
    }

    //! Create a vector with the corners of sub control volume faces
    ScvfCornerStorage getScvfCorners(unsigned int localScvfIdx) const
    {
        return ScvfCornerStorage{{p_[0]}};
    }

    //! Create the sub control volume face geometries on the boundary
    ScvfCornerStorage getBoundaryScvfCorners(const Intersection& is,
                                             const typename Intersection::Geometry& geometry,
                                             unsigned int indexInIntersection) const
    {
        return ScvfCornerStorage{{geometry.corner(0)}};
    }

    //! get scvf normal vector
    GlobalPosition normal(const ScvfCornerStorage& scvfCorners,
                          const std::vector<unsigned int>& scvIndices) const
    {
        auto normal = p_[2] - p_[1];
        normal /= normal.two_norm();
        return normal;
    }

    //! get scv volume
    Scalar scvVolume(const ScvCornerStorage& scvCorners) const
    {
        return (scvCorners[1] - scvCorners[0]).two_norm();
    }

    //! get scvf area
    Scalar scvfArea(const ScvfCornerStorage& scvfCorners) const
    {
        return 1.0;
    }

protected:
    const typename Element::Geometry& elementGeometry_; //!< Reference to the element geometry
    std::size_t corners_; // number of element corners
    std::array<GlobalPosition, maxPoints> p_; // the points needed for construction of the geometries
};

//! A class to create sub control volume and sub control volume face geometries per element
template <class GridView, class ScvType, class ScvfType>
class BoxGeometryHelper<GridView, 2, ScvType, ScvfType>
{
    using Scalar = typename GridView::ctype;
    using GlobalPosition = typename Dune::FieldVector<Scalar, GridView::dimensionworld>;
    using ScvCornerStorage = typename ScvType::Traits::CornerStorage;
    using ScvfCornerStorage = typename ScvfType::Traits::CornerStorage;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    static constexpr auto dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;

    //! the maximum number of helper points used to construct the geometries
    //! Using a statically sized point array is much faster than dynamic allocation
    static constexpr int maxPoints = 9;

public:

    BoxGeometryHelper(const typename Element::Geometry& geometry)
    : elementGeometry_(geometry)
    , corners_(geometry.corners())
    {
        const auto refElement = referenceElement(geometry);

        // the element center
        p_[0] = geometry.center();

        // vertices
        for (int i = 0; i < corners_; ++i)
            p_[i+1] = geometry.corner(i);

        // face midpoints
        for (int i = 0; i < refElement.size(1); ++i)
            p_[i+corners_+1] = geometry.global(refElement.position(i, 1));
    }

    //! Create a vector with the scv corners
    ScvCornerStorage getScvCorners(unsigned int localScvIdx) const
    {
        // proceed according to number of corners of the element
        switch (corners_)
        {
        case 3: // triangle
        {
            //! Only build the maps the first time we encounter a triangle
            static const std::uint8_t vo = 1; //!< vertex offset in point vector p
            static const std::uint8_t fo = 4; //!< face offset in point vector p
            static const std::uint8_t map[3][4] =
            {
                {vo+0, fo+0, fo+1, 0},
                {vo+1, fo+2, fo+0, 0},
                {vo+2, fo+1, fo+2, 0}
            };

            return ScvCornerStorage{ {p_[map[localScvIdx][0]],
                                      p_[map[localScvIdx][1]],
                                      p_[map[localScvIdx][2]],
                                      p_[map[localScvIdx][3]]} };
        }
        case 4: // quadrilateral
        {
            //! Only build the maps the first time we encounter a quadrilateral
            static const std::uint8_t vo = 1; //!< vertex offset in point vector p
            static const std::uint8_t fo = 5; //!< face offset in point vector p
            static const std::uint8_t map[4][4] =
            {
                {vo+0, fo+2, fo+0, 0},
                {vo+1, fo+1, fo+2, 0},
                {vo+2, fo+0, fo+3, 0},
                {vo+3, fo+3, fo+1, 0}
            };

            return ScvCornerStorage{ {p_[map[localScvIdx][0]],
                                      p_[map[localScvIdx][1]],
                                      p_[map[localScvIdx][2]],
                                      p_[map[localScvIdx][3]]} };
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Box scv geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " corners=" << corners_);
        }
    }


    //! Create a vector with the corners of sub control volume faces
    ScvfCornerStorage getScvfCorners(unsigned int localScvfIdx) const
    {
        // proceed according to number of corners
        switch (corners_)
        {
        case 3: // triangle
        {
            //! Only build the maps the first time we encounter a triangle
            static const std::uint8_t fo = 4; //!< face offset in point vector p
            static const std::uint8_t map[3][2] =
            {
                {0, fo+0},
                {fo+1, 0},
                {0, fo+2}
            };

            return ScvfCornerStorage{ {p_[map[localScvfIdx][0]],
                                       p_[map[localScvfIdx][1]]} };
        }
        case 4: // quadrilateral
        {
            //! Only build the maps the first time we encounter a quadrilateral
            static const std::uint8_t fo = 5; //!< face offset in point vector p
            static const std::uint8_t map[4][2] =
            {
                {fo+0, 0},
                {0, fo+1},
                {0, fo+2},
                {fo+3, 0}
            };

            return ScvfCornerStorage{ {p_[map[localScvfIdx][0]],
                                       p_[map[localScvfIdx][1]]} };
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Box scvf geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " corners=" << corners_);
        }
    }

    //! Create the sub control volume face geometries on the boundary
    ScvfCornerStorage getBoundaryScvfCorners(const Intersection& is,
                                             const typename Intersection::Geometry& isGeom,
                                             unsigned int indexInIntersection) const
    {
        const auto refElement = referenceElement(elementGeometry_);

        const auto vIdxLocal = refElement.subEntity(is.indexInInside(), 1, indexInIntersection, dim);
        if (indexInIntersection == 0)
            return ScvfCornerStorage({p_[vIdxLocal+1], isGeom.center()});
        else if (indexInIntersection == 1)
            return ScvfCornerStorage({isGeom.center(), p_[vIdxLocal+1]});
        else
            DUNE_THROW(Dune::InvalidStateException, "local index exceeds the number of corners of 2d intersections");
    }

    //! get scvf normal vector for dim == 2, dimworld == 3
    template <int w = dimWorld>
    typename std::enable_if<w == 3, GlobalPosition>::type
    normal(const ScvfCornerStorage& scvfCorners,
           const std::vector<unsigned int>& scvIndices) const
    {
        const auto v1 = elementGeometry_.corner(1) - elementGeometry_.corner(0);
        const auto v2 = elementGeometry_.corner(2) - elementGeometry_.corner(0);
        const auto v3 = Dumux::crossProduct(v1, v2);
        const auto t = scvfCorners[1] - scvfCorners[0];
        GlobalPosition normal = Dumux::crossProduct(v3, t);
        normal /= normal.two_norm();

        //! ensure the right direction of the normal
        const auto v = elementGeometry_.corner(scvIndices[1]) - elementGeometry_.corner(scvIndices[0]);
        const auto s = v*normal;
        if (std::signbit(s))
            normal *= -1;

        return normal;
    }

    //! get scvf normal vector for dim == 2, dimworld == 2
    template <int w = dimWorld>
    typename std::enable_if<w == 2, GlobalPosition>::type
    normal(const ScvfCornerStorage& scvfCorners,
           const std::vector<unsigned int>& scvIndices) const
    {
        //! obtain normal vector by 90° counter-clockwise rotation of t
        const auto t = scvfCorners[1] - scvfCorners[0];
        GlobalPosition normal({-t[1], t[0]});
        normal /= normal.two_norm();

        //! ensure the right direction of the normal
        const auto v = elementGeometry_.corner(scvIndices[1]) - elementGeometry_.corner(scvIndices[0]);
        const auto s = v*normal;
        if (std::signbit(s))
            normal *= -1;

        return normal;
    }

    //! get scv volume for dim == 2, dimworld == 3
    template <int w = dimWorld>
    typename std::enable_if<w == 3, Scalar>::type
    scvVolume(const ScvCornerStorage& p) const
    {
        return 0.5*Dumux::crossProduct(p[3]-p[0], p[2]-p[1]).two_norm();
    }

    //! get scv volume for dim == 2, dimworld == 2
    template <int w = dimWorld>
    typename std::enable_if<w == 2, Scalar>::type
    scvVolume(const ScvCornerStorage& p) const
    {
        //! make sure we are using positive volumes
        //! Cross product of diagonals might be negative, depending on element orientation
        using std::abs;
        return 0.5*abs(Dumux::crossProduct(p[3]-p[0], p[2]-p[1]));
    }

    //! get scvf area
    Scalar scvfArea(const ScvfCornerStorage& p) const
    {
        return (p[1]-p[0]).two_norm();
    }

protected:
    const typename Element::Geometry& elementGeometry_; //!< Reference to the element geometry
    std::size_t corners_; // number of element corners
    std::array<GlobalPosition, maxPoints> p_; // the points needed for construction of the geometries
};

//! A class to create sub control volume and sub control volume face geometries per element
template <class GridView, class ScvType, class ScvfType>
class BoxGeometryHelper<GridView, 3, ScvType, ScvfType>
{
    using Scalar = typename GridView::ctype;
    using GlobalPosition = typename Dune::FieldVector<Scalar, GridView::dimensionworld>;
    using ScvCornerStorage = typename ScvType::Traits::CornerStorage;
    using ScvfCornerStorage = typename ScvfType::Traits::CornerStorage;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    static constexpr auto dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;

    //! the maximum number of helper points used to construct the geometries
    //! Using a statically sized point array is much faster than dynamic allocation
    static constexpr int maxPoints = 27;
public:
    BoxGeometryHelper(const typename Element::Geometry& geometry)
    : elementGeometry_(geometry)
    , corners_(geometry.corners())
    {
        const auto refElement = referenceElement(geometry);

        // the element center
        p_[0] = geometry.center();

        // vertices
        for (int i = 0; i < corners_; ++i)
            p_[i+1] = geometry.corner(i);

        // edge midpoints
        for (int i = 0; i < refElement.size(dim-1); ++i)
            p_[i+corners_+1] = geometry.global(refElement.position(i, dim-1));

        // face midpoints
        for (int i = 0; i < refElement.size(1); ++i)
            p_[i+corners_+1+refElement.size(dim-1)] = geometry.global(refElement.position(i, 1));
    }

    //! Create a vector with the scv corners
    ScvCornerStorage getScvCorners(unsigned int localScvIdx) const
    {
        // proceed according to number of corners of the element
        switch (corners_)
        {
        case 4: // tetrahedron
        {
            //! Only build the maps the first time we encounter a tetrahedron
            static const std::uint8_t vo = 1; //!< vertex offset in point vector p
            static const std::uint8_t eo = 5; //!< edge offset in point vector p
            static const std::uint8_t fo = 11; //!< face offset in point vector p
            static const std::uint8_t map[4][8] =
            {
                {vo+0, eo+0, eo+1, fo+0, eo+3, fo+1, fo+2,    0},
                {vo+1, eo+2, eo+0, fo+0, eo+4, fo+3, fo+1,    0},
                {vo+2, eo+1, eo+2, fo+0, eo+5, fo+2, fo+3,    0},
                {vo+3, eo+3, eo+5, fo+2, eo+4, fo+1, fo+3,    0}
            };

            return ScvCornerStorage{ {p_[map[localScvIdx][0]],
                                      p_[map[localScvIdx][1]],
                                      p_[map[localScvIdx][2]],
                                      p_[map[localScvIdx][3]],
                                      p_[map[localScvIdx][4]],
                                      p_[map[localScvIdx][5]],
                                      p_[map[localScvIdx][6]],
                                      p_[map[localScvIdx][7]]} };
        }
        case 6: // prism
        {
            //! Only build the maps the first time we encounter a prism
            static const std::uint8_t vo = 1; //!< vertex offset in point vector p
            static const std::uint8_t eo = 7; //!< edge offset in point vector p
            static const std::uint8_t fo = 16; //!< face offset in point vector p
            static const std::uint8_t map[6][8] =
            {
                {vo+0, eo+3, eo+4, fo+3, eo+0, fo+0, fo+1,    0},
                {vo+1, eo+5, eo+3, fo+3, eo+1, fo+2, fo+0,    0},
                {vo+2, eo+4, eo+5, fo+3, eo+2, fo+1, fo+2,    0},
                {vo+3, eo+7, eo+6, fo+4, eo+0, fo+1, fo+0,    0},
                {vo+4, eo+6, eo+8, fo+4, eo+1, fo+0, fo+2,    0},
                {vo+5, eo+8, eo+7, fo+4, eo+2, fo+2, fo+1,    0}
            };

            return ScvCornerStorage{ {p_[map[localScvIdx][0]],
                                      p_[map[localScvIdx][1]],
                                      p_[map[localScvIdx][2]],
                                      p_[map[localScvIdx][3]],
                                      p_[map[localScvIdx][4]],
                                      p_[map[localScvIdx][5]],
                                      p_[map[localScvIdx][6]],
                                      p_[map[localScvIdx][7]]} };
        }
        case 8: // hexahedron
        {
            //! Only build the maps the first time we encounter a hexahedron
            static const std::uint8_t vo = 1; //!< vertex offset in point vector p
            static const std::uint8_t eo = 9; //!< edge offset in point vector p
            static const std::uint8_t fo = 21; //!< face offset in point vector p
            static const std::uint8_t map[8][8] =
            {
                {vo+0, eo+6, eo+4, fo+4, eo+0, fo+2, fo+0,    0},
                {vo+1, eo+5, eo+6, fo+4, eo+1, fo+1, fo+2,    0},
                {vo+2, eo+4, eo+7, fo+4, eo+2, fo+0, fo+3,    0},
                {vo+3, eo+7, eo+5, fo+4, eo+3, fo+3, fo+1,    0},
                {vo+4, eo+8, eo+10, fo+5, eo+0, fo+0, fo+2,    0},
                {vo+5, eo+10, eo+9, fo+5, eo+1, fo+2, fo+1,    0},
                {vo+6, eo+11, eo+8, fo+5, eo+2, fo+3, fo+0,    0},
                {vo+7, eo+9, eo+11, fo+5, eo+3, fo+1, fo+3,    0},
            };

            return ScvCornerStorage{ {p_[map[localScvIdx][0]],
                                      p_[map[localScvIdx][1]],
                                      p_[map[localScvIdx][2]],
                                      p_[map[localScvIdx][3]],
                                      p_[map[localScvIdx][4]],
                                      p_[map[localScvIdx][5]],
                                      p_[map[localScvIdx][6]],
                                      p_[map[localScvIdx][7]]} };
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Box scv geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " corners=" << corners_);
        }
    }

    //! Create a vector with the scvf corners
    ScvfCornerStorage getScvfCorners(unsigned int localScvfIdx) const
    {
        // proceed according to number of corners of the element
        switch (corners_)
        {
        case 4: // tetrahedron
        {
            //! Only build the maps the first time we encounter a tetrahedron
            static const std::uint8_t eo = 5; //!< edge offset in point vector p
            static const std::uint8_t fo = 11; //!< face offset in point vector p
            static const std::uint8_t map[6][4] =
            {
                {eo+0, fo+0, fo+1,    0},
                {fo+0, eo+1,    0, fo+2},
                {eo+2, fo+0, fo+3,    0},
                {fo+2, eo+3,    0, fo+1},
                {fo+3,    0, eo+4, fo+1},
                {eo+5, fo+2, fo+3,    0}
            };

            return ScvfCornerStorage{ {p_[map[localScvfIdx][0]],
                                       p_[map[localScvfIdx][1]],
                                       p_[map[localScvfIdx][2]],
                                       p_[map[localScvfIdx][3]]} };
        }
        case 6: // prism
        {
            //! Only build the maps the first time we encounter a prism
            static const std::uint8_t eo = 7; //!< edge offset in point vector p
            static const std::uint8_t fo = 16; //!< face offset in point vector p
            static const std::uint8_t map[9][4] =
            {
                {eo+0, fo+0, fo+1, 0},
                {eo+1, fo+2, fo+0, 0},
                {eo+2, fo+1, fo+2, 0},
                {eo+3, fo+0, fo+3, 0},
                {eo+4, fo+3, fo+1, 0},
                {eo+5, fo+2, fo+3, 0},
                {eo+6, fo+4, fo+0, 0},
                {eo+7, fo+1, fo+4, 0},
                {eo+8, fo+4, fo+2, 0}
            };

            return ScvfCornerStorage{ {p_[map[localScvfIdx][0]],
                                       p_[map[localScvfIdx][1]],
                                       p_[map[localScvfIdx][2]],
                                       p_[map[localScvfIdx][3]]} };
        }
        case 8: // hexahedron
        {
            //! Only build the maps the first time we encounter a hexahedron
            static const std::uint8_t eo = 9; //!< edge offset in point vector p
            static const std::uint8_t fo = 21; //!< face offset in point vector p
            static const std::uint8_t map[12][4] =
            {
                {fo+0, eo+0,    0, fo+2},
                {fo+1,    0, eo+1, fo+2},
                {fo+3, eo+2,    0, fo+0},
                {eo+3, fo+3, fo+1,    0},
                {fo+4, eo+4,    0, fo+0},
                {eo+5, fo+4, fo+1,    0},
                {eo+6, fo+4, fo+2,    0},
                {fo+4, eo+7,    0, fo+3},
                {   0, fo+0, fo+5, eo+8},
                {eo+9, fo+1, fo+5,    0},
                {eo+10, fo+2, fo+5,   0},
                {eo+11, fo+5, fo+3,   0}
            };

            return ScvfCornerStorage{ {p_[map[localScvfIdx][0]],
                                       p_[map[localScvfIdx][1]],
                                       p_[map[localScvfIdx][2]],
                                       p_[map[localScvfIdx][3]]} };
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Box scv geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " corners=" << corners_);
        }
    }

    //! Create the sub control volume face geometries on the boundary
    ScvfCornerStorage getBoundaryScvfCorners(const Intersection& is,
                                             const typename Intersection::Geometry& geometry,
                                             unsigned int indexInIntersection) const
    {
        const auto refElement = referenceElement(elementGeometry_);
        const auto faceRefElem = referenceElement(geometry);

        GlobalPosition pi[9];
        auto corners = geometry.corners();

        // the facet center
        pi[0] = geometry.center();

        // corners
        const auto idxInInside = is.indexInInside();
        for (int i = 0; i < corners; ++i)
        {
            const auto vIdxLocal = refElement.subEntity(idxInInside, 1, i, dim);
            pi[i+1] = elementGeometry_.corner(vIdxLocal);
        }

        // edge midpoints
        for (int i = 0; i < faceRefElem.size(1); ++i)
        {
            const auto edgeIdxLocal = refElement.subEntity(idxInInside, 1, i, dim-1);
            pi[i+corners+1] = p_[edgeIdxLocal+corners_+1];
        }

        // procees according to number of corners
        switch (corners)
        {
        case 3: // triangle
        {
            //! Only build the maps the first time we encounter a triangle
            static const std::uint8_t vo = 1; //!< vertex offset in point vector pi
            static const std::uint8_t eo = 4; //!< edge offset in point vector pi
            static const std::uint8_t map[3][4] =
            {
                {vo+0, eo+0, eo+1, 0},
                {vo+1, eo+2, eo+0, 0},
                {vo+2, eo+1, eo+2, 0}
            };

            return ScvfCornerStorage{ {pi[map[indexInIntersection][0]],
                                       pi[map[indexInIntersection][1]],
                                       pi[map[indexInIntersection][2]],
                                       pi[map[indexInIntersection][3]]} };
        }
        case 4: // quadrilateral
        {
            //! Only build the maps the first time we encounter a quadrilateral
            static const std::uint8_t vo = 1; //!< vertex offset in point vector pi
            static const std::uint8_t eo = 5; //!< edge offset in point vector pi
            static const std::uint8_t map[4][4] =
            {
                {vo+0, eo+2, eo+0, 0},
                {vo+1, eo+1, eo+2, 0},
                {vo+2, eo+0, eo+3, 0},
                {vo+3, eo+3, eo+1, 0}
            };

            return ScvfCornerStorage{ {pi[map[indexInIntersection][0]],
                                       pi[map[indexInIntersection][1]],
                                       pi[map[indexInIntersection][2]],
                                       pi[map[indexInIntersection][3]]} };
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Box scvf boundary geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " corners=" << corners);
        }
    }

    //! get scvf normal vector
    GlobalPosition normal(const ScvfCornerStorage& p,
                          const std::vector<unsigned int>& scvIndices) const
    {
        auto normal = Dumux::crossProduct(p[1]-p[0], p[2]-p[0]);
        normal /= normal.two_norm();

        const auto v = elementGeometry_.corner(scvIndices[1]) - elementGeometry_.corner(scvIndices[0]);
        const auto s = v*normal;
        if (std::signbit(s))
            normal *= -1;

        return normal;
    }

    //! get scv volume
    Scalar scvVolume(const ScvCornerStorage& p) const
    {
        // after Grandy 1997, Efficient computation of volume of hexahedron
        const auto v = p[7]-p[0];
        return 1.0/6.0 * ( Dumux::tripleProduct(v, p[1]-p[0], p[3]-p[5])
                         + Dumux::tripleProduct(v, p[4]-p[0], p[5]-p[6])
                         + Dumux::tripleProduct(v, p[2]-p[0], p[6]-p[3]));
    }

    //! get scvf area
    Scalar scvfArea(const ScvfCornerStorage& p) const
    {
        // after Wolfram alpha quadrilateral area
        return 0.5*Dumux::crossProduct(p[3]-p[0], p[2]-p[1]).two_norm();
    }

protected:
    const typename Element::Geometry& elementGeometry_; //!< Reference to the element geometry
    std::size_t corners_; // number of element corners
    std::array<GlobalPosition, maxPoints> p_; // the points needed for construction of the scv/scvf geometries
};

} // end namespace Dumux

#endif
