// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Helper class constructing the dual grid finite volume geometries
 *        for the box discretizazion method
 */
#ifndef DUMUX_DISCRETIZATION_BOX_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_BOX_GEOMETRY_HELPER_HH

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/implicit/box/properties.hh>
#include <dumux/common/math.hh>

namespace Dumux
{

//! Create sub control volumes and sub control volume face geometries
template<class GridView, int dim>
class BoxGeometryHelper;

//! A class to create sub control volume and sub control volume face geometries per element
template <class GridView>
class BoxGeometryHelper<GridView, 1>
{
private:
    using Scalar = typename GridView::ctype;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using ScvGeometry = Dune::CachedMultiLinearGeometry<Scalar, dim, dimWorld>;
    using ScvfGeometry = Dune::CachedMultiLinearGeometry<Scalar, dim-1, dimWorld>;

    using GlobalPosition = typename ScvGeometry::GlobalCoordinate;
    using CornerList = std::vector<GlobalPosition>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;

public:
    using ScvGeometryVector = std::vector<ScvGeometry>;
    using ScvfGeometryVector = std::vector<ScvfGeometry>;

    //! get sub control volume geometries from element
    void getScvAndScvfGeometries(const typename Element::Geometry& geometry,
                                 ScvGeometryVector& scvGeometries,
                                 ScvfGeometryVector& scvfGeometries)
    {
        scvGeometries.reserve(2);
        scvfGeometries.reserve(1);

        // the sub control volumes
        scvGeometries.emplace_back(Dune::GeometryType(1), std::vector<GlobalPosition>({geometry.corner(0), geometry.center()}));
        scvGeometries.emplace_back(Dune::GeometryType(1), std::vector<GlobalPosition>({geometry.center(), geometry.corner(1)}));

        // the sub control volume faces
        scvfGeometries.emplace_back(Dune::GeometryType(0), std::vector<GlobalPosition>({geometry.center()}));
    }

    //! get sub control volume geometries from element of dimension 1
    ScvfGeometryVector getBoundaryScvfGeometries(const typename Intersection::Geometry& geometry)
    {
        return {ScvfGeometry(Dune::GeometryType(0), std::vector<GlobalPosition>({geometry.center()}))};
    }

    //! get scvf normal vector
    static GlobalPosition normal(const typename Element::Geometry& geometry,
                                 const ScvfGeometry& scvfGeometry)
    {
        GlobalPosition normal = geometry.corner(1) - geometry.corner(0);
        normal /= normal.two_norm();
        return normal;
    }
};

//! A class to create sub control volume and sub control volume face geometries per element
template <class GridView>
class BoxGeometryHelper<GridView, 2>
{
    using Scalar = typename GridView::ctype;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using ScvGeometry = Dune::CachedMultiLinearGeometry<Scalar, dim, dimWorld>;
    using ScvfGeometry = Dune::CachedMultiLinearGeometry<Scalar, dim-1, dimWorld>;

    using GlobalPosition = typename ScvGeometry::GlobalCoordinate;
    using PointVector = std::vector<GlobalPosition>;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;

    //! the maximum number of helper points used to construct the geometries
    //! Using a statically sized point array is much faster than dynamic allocation
    static constexpr int maxPoints = 9;

public:

    BoxGeometryHelper(const typename Element::Geometry& geometry)
    : elementGeometry_(geometry), corners_(geometry.corners())
    {
        // extract the corners of the sub control volumes
        const auto& referenceElement = ReferenceElements::general(geometry.type());

        // the element center
        p[0] = geometry.center();

        // vertices
        for (int i = 0; i < corners_; ++i)
            p[i+1] = geometry.corner(i);

        // face midpoints
        for (int i = 0; i < referenceElement.size(1); ++i)
            p[i+corners_+1] = geometry.global(referenceElement.position(i, 1));
    }

    //! Create a vector with the scv corners
    PointVector getScvCorners(unsigned int localScvIdx) const
    {
        // proceed according to number of corners of the element
        switch (corners_)
        {
        case 3: // triangle
        {
            //! Only build the maps the first time we encounter a triangle
            static const std::uint8_t vo = 1; //! vertex offset in point vector p
            static const std::uint8_t fo = 4; //! face offset in point vector p
            static const std::uint8_t map[3][4] =
            {
                {vo+0, fo+0, fo+1, 0},
                {vo+1, fo+2, fo+0, 0},
                {vo+2, fo+1, fo+2, 0}
            };

            return PointVector( {p[map[localScvIdx][0]],
                                 p[map[localScvIdx][1]],
                                 p[map[localScvIdx][2]],
                                 p[map[localScvIdx][3]]} );
        }
        case 4: // quadrilateral
        {
            //! Only build the maps the first time we encounter a quadrilateral
            static const std::uint8_t vo = 1; //! vertex offset in point vector p
            static const std::uint8_t fo = 5; //! face offset in point vector p
            static const std::uint8_t map[4][4] =
            {
                {vo+0, fo+2, fo+0, 0},
                {vo+1, fo+1, fo+2, 0},
                {vo+2, fo+0, fo+3, 0},
                {vo+3, fo+3, fo+1, 0}
            };

            return PointVector( {p[map[localScvIdx][0]],
                                 p[map[localScvIdx][1]],
                                 p[map[localScvIdx][2]],
                                 p[map[localScvIdx][3]]} );
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Box scv geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " corners=" << corners_);
        }
    }


    //! Create a vector with the corners of sub control volume faces
    PointVector getScvfCorners(unsigned int localScvfIdx) const
    {
        // proceed according to number of corners
        switch (corners_)
        {
        case 3: // triangle
        {
            //! Only build the maps the first time we encounter a triangle
            static const std::uint8_t fo = 4; //! face offset in point vector p
            static const std::uint8_t map[3][2] =
            {
                {0, fo+0},
                {fo+1, 0},
                {0, fo+2}
            };

            return PointVector( {p[map[localScvfIdx][0]],
                                 p[map[localScvfIdx][1]]} );
        }
        case 4: // quadrilateral
        {
            //! Only build the maps the first time we encounter a quadrilateral
            static const std::uint8_t fo = 5; //! face offset in point vector p
            static const std::uint8_t map[4][2] =
            {
                {fo+0, 0},
                {0, fo+1},
                {0, fo+2},
                {fo+3, 0}
            };

            return PointVector( {p[map[localScvfIdx][0]],
                                 p[map[localScvfIdx][1]]} );
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Box scvf geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " corners=" << corners_);
        }
    }

    //! Create the sub control volume face geometries on the boundary
    PointVector getBoundaryScvfCorners(const typename Intersection::Geometry& geometry,
                                       unsigned int indexInIntersection) const
    {
        if (indexInIntersection == 0)
            return PointVector({geometry.corner(0), geometry.center()});
        else if (indexInIntersection == 1)
            return PointVector({geometry.center(), geometry.corner(1)});
        else
            DUNE_THROW(Dune::InvalidStateException, "local index exceeds the number of corners of 2d intersections");
    }

    //! get scvf normal vector for dim == 2, dimworld == 3
    template <int w = dimWorld>
    typename std::enable_if<w == 3, GlobalPosition>::type
    normal(const PointVector& scvfCorners,
           const std::vector<unsigned int>& scvIndices) const
    {
        const auto v1 = elementGeometry_.corner(1) - elementGeometry_.corner(0);
        const auto v2 = elementGeometry_.corner(2) - elementGeometry_.corner(0);
        const auto v3 = Dumux::crossProduct(v1, v2);
        const auto t = scvfCorners[1] - scvfCorners[0];
        GlobalPosition normal = Dumux::crossProduct(v3, t);
        normal /= normal.two_norm();

        // TODO can this be done easier?, e.g. always ensure the right direction?
        const auto v = elementGeometry_.corner(scvIndices[1]) - elementGeometry_.corner(scvIndices[0]);
        const auto s = v*normal;
        if (std::signbit(s))
            normal *= -1;

        return normal;
    }

    //! get scvf normal vector for dim == 2, dimworld == 2
    template <int w = dimWorld>
    typename std::enable_if<w == 2, GlobalPosition>::type
    normal(const PointVector& scvfCorners,
           const std::vector<unsigned int>& scvIndices) const
    {
        const auto t = scvfCorners[1] - scvfCorners[0];
        GlobalPosition normal({-t[1], t[0]});
        normal /= normal.two_norm();
        return normal;
    }

    //! get scv volume for dim == 2, dimworld == 3
    template <int w = dimWorld>
    typename std::enable_if<w == 3, Scalar>::type
    scvVolume(const PointVector& scvCorners) const
    {
        const auto v1 = scvCorners[1] - scvCorners[0];
        const auto v2 = scvCorners[2] - scvCorners[0];
        return Dumux::crossProduct(v1, v2).two_norm();
    }

    //! get scv volume for dim == 2, dimworld == 2
    template <int w = dimWorld>
    typename std::enable_if<w == 2, Scalar>::type
    scvVolume(const PointVector& scvCorners) const
    {
        const auto v1 = scvCorners[1] - scvCorners[0];
        const auto v2 = scvCorners[2] - scvCorners[0];
        return Dumux::crossProduct(v1, v2);
    }

    //! get scvf area
    Scalar scvfArea(const PointVector& scvfCorners) const
    {
        return (scvfCorners[1] - scvfCorners[0]).two_norm();
    }

private:
    const typename Element::Geometry& elementGeometry_; //! Reference to the element geometry
    GlobalPosition p[maxPoints]; // the points needed for construction of the geometries
    std::size_t corners_; // number of element corners
};

//! A class to create sub control volume and sub control volume face geometries per element
template <class GridView>
class BoxGeometryHelper<GridView, 3>
{
private:
    using Scalar = typename GridView::ctype;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using ScvGeometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld>;
    using ScvfGeometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld>;

    using GlobalPosition = typename ScvGeometry::GlobalCoordinate;
    using CornerList = std::vector<GlobalPosition>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;
    using FaceReferenceElements = typename Dune::ReferenceElements<Scalar, dim-1>;

public:
    using ScvGeometryVector = std::vector<ScvGeometry>;
    using ScvfGeometryVector = std::vector<ScvfGeometry>;

    //! get sub control volume geometries from element
    void getScvAndScvfGeometries(const typename Element::Geometry& geometry,
                            ScvGeometryVector& scvGeometries,
                            ScvfGeometryVector& scvfGeometries)
    {
        // sub control volume geometries in 3D are always hexahedrons
        Dune::GeometryType scvGeometryType; scvGeometryType.makeHexahedron();
        // sub control volume face geometries in 3D are always quadrilaterals
        Dune::GeometryType scvfGeometryType; scvfGeometryType.makeQuadrilateral();

        // extract the corners of the sub control volumes
        const auto& referenceElement = ReferenceElements::general(geometry.type());

        // vertices
        CornerList v;
        for (int i = 0; i < geometry.corners(); ++i)
            v.emplace_back(geometry.corner(i));

        // edge midpoints
        CornerList e;
        for (int i = 0; i < referenceElement.size(dim-1); ++i)
            e.emplace_back(geometry.global(referenceElement.position(i, dim-1)));

        // face midpoints
        CornerList f;
        for (int i = 0; i < referenceElement.size(1); ++i)
            f.emplace_back(geometry.global(referenceElement.position(i, 1)));

        auto c = geometry.center();

        // procees according to number of corners
        // \todo prisms (corners == 6) and pyramids (corners == 5)
        switch (geometry.corners())
        {
        case 4: // tetrahedron
        {
            scvGeometries.reserve(4);
            scvfGeometries.reserve(6);

            // sub control volumes
            ScvGeometryVector scvGeometries(4);
            scvGeometries.emplace_back(scvGeometryType, std::vector<GlobalPosition>({v[0], e[0], e[1], f[0], e[3], f[1], f[2], c}));
            scvGeometries.emplace_back(scvGeometryType, std::vector<GlobalPosition>({v[1], e[2], e[0], f[0], f[3], e[4], c, f[1]}));
            scvGeometries.emplace_back(scvGeometryType, std::vector<GlobalPosition>({v[2], e[1], e[2], f[0], e[5], f[2], f[3], c}));
            scvGeometries.emplace_back(scvGeometryType, std::vector<GlobalPosition>({v[3], e[3], e[5], f[2], e[4], f[1], f[3], c}));

            // sub control volume faces
            ScvfGeometryVector scvfGeometries(6);
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({e[0], f[0], f[1], c}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({f[0], e[1], c, f[2]}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({e[2], f[0], f[3], c}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({f[2], e[3], c, f[1]}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({f[3], c, e[4], f[1]}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({e[5], f[2], f[3], c}));

            break;
        }
        case 8: // hexahedron
        {
            scvGeometries.reserve(8);
            scvfGeometries.reserve(12);

            // sub control volumes
            ScvGeometryVector scvGeometries(8);
            scvGeometries.emplace_back(scvGeometryType, std::vector<GlobalPosition>({v[0], e[6], e[4], f[4], e[0], f[2], f[0], c}));
            scvGeometries.emplace_back(scvGeometryType, std::vector<GlobalPosition>({e[6], v[1], f[4], e[5], f[2], e[1], c, f[1]}));
            scvGeometries.emplace_back(scvGeometryType, std::vector<GlobalPosition>({e[4], f[4], v[2], e[7], f[0], c, e[2], f[3]}));
            scvGeometries.emplace_back(scvGeometryType, std::vector<GlobalPosition>({f[4], e[5], e[7], v[3], c, f[1], f[3], e[3]}));
            scvGeometries.emplace_back(scvGeometryType, std::vector<GlobalPosition>({e[0], f[2], f[0], c, v[4], e[10], e[8], f[5]}));
            scvGeometries.emplace_back(scvGeometryType, std::vector<GlobalPosition>({f[2], e[1], c, f[1], e[10], v[5], f[5], e[9]}));
            scvGeometries.emplace_back(scvGeometryType, std::vector<GlobalPosition>({f[0], c, e[2], f[3], e[8], f[5], v[6], e[11]}));
            scvGeometries.emplace_back(scvGeometryType, std::vector<GlobalPosition>({c, f[1], f[3], e[3], f[5], e[9], e[11], v[7]}));

            // sub control volume faces
            ScvfGeometryVector scvfGeometries(12);
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({f[0], e[0], c, f[2]}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({f[1], c, e[1], f[2]}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({f[3], e[2], c, f[0]}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({e[3], f[3], f[1], c}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({f[4], e[4], c, f[0]}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({e[5], f[4], f[1], c}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({e[6], f[4], f[2], c}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({f[4], e[7], c, f[3]}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({c, f[0], f[5], e[8]}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({e[9], f[1], f[5], c}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({e[10], f[2], f[5], c}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({e[11], f[5], f[3], c}));

            break;
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Box scv geometries for " << geometry.type());
        }
    }

    //! get sub control volume geometries from element
    ScvfGeometryVector getBoundaryScvfGeometries(const typename Intersection::Geometry& geometry)
    {
        ScvfGeometryVector scvfGeometries;

        // sub control volume face geometries in 3D are always quadrilaterals
        Dune::GeometryType scvfGeometryType; scvfGeometryType.makeQuadrilateral();

        // extract the corners of the sub control volumes
        const auto& referenceElement = FaceReferenceElements::general(geometry.type());

        // vertices
        CornerList v;
        for (int i = 0; i < geometry.corners(); ++i)
            v.emplace_back(geometry.corner(i));

        // edge midpoints
        CornerList e;
        for (int i = 0; i < referenceElement.size(1); ++i)
            e.emplace_back(geometry.global(referenceElement.position(i, 1)));

        // face midpoint
        auto c = geometry.center();

        // procees according to number of corners
        switch (geometry.corners())
        {
        case 3: // triangle
        {
            scvfGeometries.reserve(3);
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({v[0], e[0], e[1], c}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({v[1], e[2], e[0], c}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({v[2], e[1], e[2], c}));

            return scvfGeometries;
        }
        case 4: // quadrilateral
        {
            scvfGeometries.reserve(4);
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({v[0], e[2], e[0], c}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({v[1], e[1], e[2], c}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({v[2], e[0], e[3], c}));
            scvfGeometries.emplace_back(scvfGeometryType, std::vector<GlobalPosition>({v[3], e[3], e[1], c}));

            return scvfGeometries;
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Box scvf boundary geometries for " << geometry.type());
        }
    }

    //! get scvf normal vector
    static GlobalPosition normal(const typename Element::Geometry& geometry,
                                 const ScvfGeometry& scvfGeometry)
    {
        const auto t1 = scvfGeometry.corner(1) - scvfGeometry.corner(0);
        const auto t2 = scvfGeometry.corner(2) - scvfGeometry.corner(0);
        GlobalPosition normal = Dumux::crossProduct(t1, t2);
        normal /= normal.two_norm();
        return normal;
    }
};

} // end namespace Dumux

#endif
