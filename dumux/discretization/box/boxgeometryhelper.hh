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
private:
    using Scalar = typename GridView::ctype;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using ScvGeometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld>;
    using ScvfGeometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld>;

    using GlobalPosition = typename ScvGeometry::GlobalCoordinate;
    using PointVector = std::vector<GlobalPosition>;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;

public:

    BoxGeometryHelper(const typename Element::Geometry& geometry)
    : vertexOffset(1),
      faceOffset(1 + geometry.corners())
    {
        // extract the corners of the sub control volumes
        const auto& referenceElement = ReferenceElements::general(geometry.type());

        // the element center
        p.emplace_back(geometry.center());

        // vertices
        for (int i = 0; i < geometry.corners(); ++i)
            p.emplace_back(geometry.corner(i));

        // face midpoints
        for (int i = 0; i < referenceElement.size(1); ++i)
            p.emplace_back(geometry.global(referenceElement.position(i, 1)));

        corners = geometry.corners();
    }

    //! Create a vector with the scv corners
    PointVector getScvCorners(unsigned int localScvIdx)
    {
        // proceed according to number of corners of the element
        switch (corners)
        {
        case 3: // triangle
        {
            static const unsigned int triangleScvMap[3][4] =
            {
                {vertexOffset+0, faceOffset+0, faceOffset+1, 0},
                {vertexOffset+1, faceOffset+2, faceOffset+0, 0},
                {vertexOffset+2, faceOffset+1, faceOffset+2, 0}
            };

            return PointVector( {p[triangleScvMap[localScvIdx][0]],
                                 p[triangleScvMap[localScvIdx][1]],
                                 p[triangleScvMap[localScvIdx][2]],
                                 p[triangleScvMap[localScvIdx][3]]} );
        }
        case 4: // quadrilateral
        {
            static const unsigned int quadScvMap[4][4] =
            {
                {vertexOffset+0, faceOffset+2, faceOffset+0, 0},
                {vertexOffset+1, faceOffset+1, faceOffset+2, 0},
                {vertexOffset+2, faceOffset+0, faceOffset+3, 0},
                {vertexOffset+3, faceOffset+3, faceOffset+1, 0}
            };

            return PointVector( {p[quadScvMap[localScvIdx][0]],
                                 p[quadScvMap[localScvIdx][1]],
                                 p[quadScvMap[localScvIdx][2]],
                                 p[quadScvMap[localScvIdx][3]]} );
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Box scv geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " corners=" << corners);
        }
    }


    //! Create a vector with the corners of sub control volume faces
    PointVector getScvfCorners(unsigned int localScvfIdx)
    {
        // proceed according to number of corners
        switch (corners)
        {
        case 3: // triangle
        {
            static const unsigned int triangleScvfMap[3][2] =
            {
                {0, faceOffset+0},
                {faceOffset+1, 0},
                {0, faceOffset+2}
            };

            return PointVector( {p[triangleScvfMap[localScvfIdx][0]],
                                 p[triangleScvfMap[localScvfIdx][1]]} );
        }
        case 4: // quadrilateral
        {
            static const unsigned int quadScvfMap[4][2] =
            {
                {faceOffset+0, 0},
                {0, faceOffset+1},
                {0, faceOffset+2},
                {faceOffset+3, 0}
            };

            return PointVector( {p[quadScvfMap[localScvfIdx][0]],
                                 p[quadScvfMap[localScvfIdx][1]]} );
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Box scvf geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " corners=" << corners);
        }
    }

    //! Create the sub control volume face geometries on the boundary
    static PointVector getBoundaryScvfCorners(const typename Intersection::Geometry& geometry,
                                                            unsigned int idxOnIntersection)
    {
        if (idxOnIntersection == 0)
            return PointVector({geometry.corner(0), geometry.center()});
        else if (idxOnIntersection == 1)
            return PointVector({geometry.center(), geometry.corner(1)});
        else
            DUNE_THROW(Dune::InvalidStateException, "local index exceeds the number of corners of 2d intersections");
    }

    //! get scvf normal vector for dim == 2, dimworld == 3
    template <int w = dimWorld>
    static typename std::enable_if<w == 3, GlobalPosition>::type
    normal(const typename Element::Geometry& geometry,
           const PointVector& scvfCorners)
    {
        const auto v1 = geometry.corner(1) - geometry.corner(0);
        const auto v2 = geometry.corner(2) - geometry.corner(0);
        const auto v3 = Dumux::crossProduct(v1, v2);
        const auto t = scvfCorners[1] - scvfCorners[0];
        GlobalPosition normal = Dumux::crossProduct(v3, t);
        normal /= normal.two_norm();
        return normal;
    }

    //! get scvf normal vector for dim == 2, dimworld == 2
    template <int w = dimWorld>
    static typename std::enable_if<w == 2, GlobalPosition>::type
    normal(const typename Element::Geometry& geometry,
           const PointVector& scvfCorners)
    {
        const auto t = scvfCorners[1] - scvfCorners[0];
        GlobalPosition normal({-t[1], t[0]});
        normal /= normal.two_norm();
        return normal;
    }

    //! get scv volume for dim == 2, dimworld == 3
    template <int w = dimWorld>
    static typename std::enable_if<w == 3, Scalar>::type
    volume(const PointVector& scvCorners)
    {
        const auto v1 = scvCorners[1] - scvCorners[0];
        const auto v2 = scvCorners[2] - scvCorners[0];
        return Dumux::crossProduct(v1, v2).two_norm();
    }

    //! get scv volume for dim == 2, dimworld == 2
    template <int w = dimWorld>
    static typename std::enable_if<w == 2, Scalar>::type
    volume(const PointVector& scvCorners)
    {
        const auto v1 = scvCorners[1] - scvCorners[0];
        const auto v2 = scvCorners[2] - scvCorners[0];
        return Dumux::crossProduct(v1, v2);
    }

    //! get scvf area
    static Scalar area(const PointVector& scvfCorners)
    {
        return (scvfCorners[1] - scvfCorners[0]).two_norm();
    }
private:
    std::vector<GlobalPosition> p; // the points needed for construction of the geometries
    std::size_t corners; // number of element corners
    const unsigned int vertexOffset; // offset for vertex positions in point vector
    const unsigned int faceOffset; // offset for face center positions in point vector
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
