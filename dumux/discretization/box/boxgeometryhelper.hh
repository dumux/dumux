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

#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/math.hh>

namespace Dumux {

//! Traits for an efficient corner storage for box method sub control volumes
template <class ct>
struct BoxMLGeometryTraits : public Dune::MultiLinearGeometryTraits<ct>
{
    // we use static vectors to store the corners as we know
    // the number of corners in advance (2^(mydim) corners (1<<(mydim))
    template< int mydim, int cdim >
    struct CornerStorage
    {
        using Type = std::array< Dune::FieldVector< ct, cdim >, (1<<(mydim)) >;
    };

    // we know all scvfs will have the same geometry type
    template< int mydim >
    struct hasSingleGeometryType
    {
        static const bool v = true;
        static const unsigned int topologyId = Dune::GeometryTypes::cube(mydim).id();
    };
};

namespace Detail {

template<Dune::GeometryType::Id gt>
struct ScvCorners;

template<>
struct ScvCorners<Dune::GeometryTypes::line>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 2>, 2> keys = {{
        { Key{0, 1}, Key{0, 0} },
        { Key{1, 1}, Key{0, 0} }
    }};
};

template<>
struct ScvCorners<Dune::GeometryTypes::triangle>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 4>, 3> keys = {{
        { Key{0, 2}, Key{0, 1}, Key{1, 1}, Key{0, 0} },
        { Key{1, 2}, Key{2, 1}, Key{0, 1}, Key{0, 0} },
        { Key{2, 2}, Key{1, 1}, Key{2, 1}, Key{0, 0} }
    }};
};

template<>
struct ScvCorners<Dune::GeometryTypes::quadrilateral>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 4>, 4> keys = {{
        { Key{0, 2}, Key{2, 1}, Key{0, 1}, Key{0, 0} },
        { Key{1, 2}, Key{1, 1}, Key{2, 1}, Key{0, 0} },
        { Key{2, 2}, Key{0, 1}, Key{3, 1}, Key{0, 0} },
        { Key{3, 2}, Key{3, 1}, Key{1, 1}, Key{0, 0} }
    }};
};

template<>
struct ScvCorners<Dune::GeometryTypes::tetrahedron>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 8>, 4> keys = {{
        { Key{0, 3}, Key{0, 2}, Key{1, 2}, Key{0, 1}, Key{3, 2}, Key{1, 1}, Key{2, 1}, Key{0, 0} },
        { Key{1, 3}, Key{2, 2}, Key{0, 2}, Key{0, 1}, Key{4, 2}, Key{3, 1}, Key{1, 1}, Key{0, 0} },
        { Key{2, 3}, Key{1, 2}, Key{2, 2}, Key{0, 1}, Key{5, 2}, Key{2, 1}, Key{3, 1}, Key{0, 0} },
        { Key{3, 3}, Key{3, 2}, Key{5, 2}, Key{2, 1}, Key{4, 2}, Key{1, 1}, Key{3, 1}, Key{0, 0} }
    }};
};

template<>
struct ScvCorners<Dune::GeometryTypes::prism>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 8>, 6> keys = {{
        { Key{0, 3}, Key{3, 2}, Key{4, 2}, Key{3, 1}, Key{0, 2}, Key{0, 1}, Key{1, 1}, Key{0, 0} },
        { Key{1, 3}, Key{5, 2}, Key{3, 2}, Key{3, 1}, Key{1, 2}, Key{2, 1}, Key{0, 1}, Key{0, 0} },
        { Key{2, 3}, Key{4, 2}, Key{5, 2}, Key{3, 1}, Key{2, 2}, Key{1, 1}, Key{2, 1}, Key{0, 0} },
        { Key{3, 3}, Key{7, 2}, Key{6, 2}, Key{4, 1}, Key{0, 2}, Key{1, 1}, Key{0, 1}, Key{0, 0} },
        { Key{4, 3}, Key{6, 2}, Key{8, 2}, Key{4, 1}, Key{1, 2}, Key{0, 1}, Key{2, 1}, Key{0, 0} },
        { Key{5, 3}, Key{8, 2}, Key{7, 2}, Key{4, 1}, Key{2, 2}, Key{2, 1}, Key{1, 1}, Key{0, 0} }
    }};
};

template<>
struct ScvCorners<Dune::GeometryTypes::hexahedron>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 8>, 8> keys = {{
        { Key{0, 3}, Key{6, 2}, Key{4, 2}, Key{4, 1}, Key{0, 2}, Key{2, 1}, Key{0, 1}, Key{0, 0} },
        { Key{1, 3}, Key{5, 2}, Key{6, 2}, Key{4, 1}, Key{1, 2}, Key{1, 1}, Key{2, 1}, Key{0, 0} },
        { Key{2, 3}, Key{4, 2}, Key{7, 2}, Key{4, 1}, Key{2, 2}, Key{0, 1}, Key{3, 1}, Key{0, 0} },
        { Key{3, 3}, Key{7, 2}, Key{5, 2}, Key{4, 1}, Key{3, 2}, Key{3, 1}, Key{1, 1}, Key{0, 0} },
        { Key{4, 3}, Key{8, 2}, Key{10, 2}, Key{5, 1}, Key{0, 2}, Key{0, 1}, Key{2, 1}, Key{0, 0} },
        { Key{5, 3}, Key{10, 2}, Key{9, 2}, Key{5, 1}, Key{1, 2}, Key{2, 1}, Key{1, 1}, Key{0, 0} },
        { Key{6, 3}, Key{11, 2}, Key{8, 2}, Key{5, 1}, Key{2, 2}, Key{3, 1}, Key{0, 1}, Key{0, 0} },
        { Key{7, 3}, Key{9, 2}, Key{11, 2}, Key{5, 1}, Key{3, 2}, Key{1, 1}, Key{3, 1}, Key{0, 0} }
    }};
};

template<Dune::GeometryType::Id gt>
struct ScvfCorners;

template<>
struct ScvfCorners<Dune::GeometryTypes::line>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 1>, 1> keys = {{
        { Key{0, 0} }
    }};
};

template<>
struct ScvfCorners<Dune::GeometryTypes::triangle>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 2>, 3> keys = {{
        { Key{0, 0}, Key{0, 1} },
        { Key{1, 1}, Key{0, 0} },
        { Key{0, 0}, Key{2, 1} }
    }};
};

template<>
struct ScvfCorners<Dune::GeometryTypes::quadrilateral>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 2>, 4> keys = {{
        { Key{0, 1}, Key{0, 0} },
        { Key{0, 0}, Key{1, 1} },
        { Key{0, 0}, Key{2, 1} },
        { Key{3, 1}, Key{0, 0} }
    }};
};

template<>
struct ScvfCorners<Dune::GeometryTypes::tetrahedron>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 4>, 6> keys = {{
        { Key{0, 2}, Key{0, 1}, Key{1, 1}, Key{0, 0} },
        { Key{0, 1}, Key{1, 2}, Key{0, 0}, Key{2, 1} },
        { Key{2, 2}, Key{0, 1}, Key{3, 1}, Key{0, 0} },
        { Key{2, 1}, Key{3, 2}, Key{0, 0}, Key{1, 1} },
        { Key{3, 1}, Key{0, 0}, Key{4, 2}, Key{1, 1} },
        { Key{5, 2}, Key{2, 1}, Key{3, 1}, Key{0, 0} }
    }};
};

template<>
struct ScvfCorners<Dune::GeometryTypes::prism>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 4>, 9> keys = {{
        { Key{0, 2}, Key{0, 1}, Key{1, 1}, Key{0, 0} },
        { Key{1, 2}, Key{2, 1}, Key{0, 1}, Key{0, 0} },
        { Key{2, 2}, Key{1, 1}, Key{2, 1}, Key{0, 0} },
        { Key{3, 2}, Key{0, 1}, Key{3, 1}, Key{0, 0} },
        { Key{4, 2}, Key{3, 1}, Key{1, 1}, Key{0, 0} },
        { Key{5, 2}, Key{2, 1}, Key{3, 1}, Key{0, 0} },
        { Key{6, 2}, Key{4, 1}, Key{0, 1}, Key{0, 0} },
        { Key{7, 2}, Key{1, 1}, Key{4, 1}, Key{0, 0} },
        { Key{8, 2}, Key{4, 1}, Key{2, 1}, Key{0, 0} }
    }};
};

template<>
struct ScvfCorners<Dune::GeometryTypes::hexahedron>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 4>, 12> keys = {{
        { Key{0, 1}, Key{0, 2}, Key{0, 0}, Key{2, 1} },
        { Key{1, 1}, Key{0, 0}, Key{1, 2}, Key{2, 1} },
        { Key{3, 1}, Key{2, 2}, Key{0, 0}, Key{0, 1} },
        { Key{3, 2}, Key{3, 1}, Key{1, 1}, Key{0, 0} },
        { Key{4, 1}, Key{4, 2}, Key{0, 0}, Key{0, 1} },
        { Key{5, 2}, Key{4, 1}, Key{1, 1}, Key{0, 0} },
        { Key{6, 2}, Key{4, 1}, Key{2, 1}, Key{0, 0} },
        { Key{4, 1}, Key{7, 2}, Key{0, 0}, Key{3, 1} },
        { Key{0, 0}, Key{0, 1}, Key{5, 1}, Key{8, 2} },
        { Key{9, 2}, Key{1, 1}, Key{5, 1}, Key{0, 0} },
        { Key{10, 2}, Key{2, 1}, Key{5, 1}, Key{0, 0} },
        { Key{11, 2}, Key{5, 1}, Key{3, 1}, Key{0, 0} }
    }};
};

// convert key array to global corner storage
template<class S, class Geo, class KeyArray, std::size_t... I>
S keyToCornerStorageImpl(const Geo& geo, const KeyArray& key, std::index_sequence<I...>)
{
    using Dune::referenceElement;
    const auto ref = referenceElement(geo);
    return { geo.global(ref.position(key[I].first, key[I].second))... };
}

// convert key array to global corner storage
template<class S, class Geo, class T, std::size_t N, class Indices = std::make_index_sequence<N>>
S keyToCornerStorage(const Geo& geo, const std::array<T, N>& key)
{
    return keyToCornerStorageImpl<S>(geo, key, Indices{});
}

} // end namespace Detail

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
    using ScvGeometry = typename ScvType::Traits::Geometry;
    using ScvfGeometry = typename ScvfType::Traits::Geometry;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
public:

    BoxGeometryHelper(const typename Element::Geometry& geometry)
    : geo_(geometry)
    {}

    //! Create a vector with the scv corners
    ScvCornerStorage getScvCorners(unsigned int localScvIdx) const
    {
        using Corners = Detail::ScvCorners<Dune::GeometryTypes::line>;
        return Detail::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localScvIdx]);
    }

    ScvGeometry scvGeometry(unsigned int localScvIdx) const
    {
        return { Dune::GeometryTypes::line, getScvCorners(localScvIdx) };
    }

    //! Create a vector with the corners of sub control volume faces
    ScvfCornerStorage getScvfCorners(unsigned int localScvfIdx) const
    {
        using Corners = Detail::ScvfCorners<Dune::GeometryTypes::line>;
        return Detail::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localScvfIdx]);
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
        auto normal = geo_.corner(1) - geo_.corner(0);
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

    //! the wrapped element geometry
    const typename Element::Geometry& elementGeometry() const
    { return geo_; }

private:
    const typename Element::Geometry& geo_; //!< Reference to the element geometry
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
public:

    BoxGeometryHelper(const typename Element::Geometry& geometry)
    : geo_(geometry)
    {}

    //! Create a vector with the scv corners
    ScvCornerStorage getScvCorners(unsigned int localScvIdx) const
    {
        // proceed according to number of corners of the element
        const auto type = geo_.type();
        if (type == Dune::GeometryTypes::triangle)
        {
            using Corners = Detail::ScvCorners<Dune::GeometryTypes::triangle>;
            return Detail::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localScvIdx]);
        }
        else if (type == Dune::GeometryTypes::quadrilateral)
        {
            using Corners = Detail::ScvCorners<Dune::GeometryTypes::quadrilateral>;
            return Detail::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localScvIdx]);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Box scv geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " type=" << type);
    }

    //! Create a vector with the corners of sub control volume faces
    ScvfCornerStorage getScvfCorners(unsigned int localScvfIdx) const
    {
        // proceed according to number of corners
        const auto type = geo_.type();
        if (type == Dune::GeometryTypes::triangle)
        {
            using Corners = Detail::ScvfCorners<Dune::GeometryTypes::triangle>;
            return Detail::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localScvfIdx]);
        }
        else if (type == Dune::GeometryTypes::quadrilateral)
        {
            using Corners = Detail::ScvfCorners<Dune::GeometryTypes::quadrilateral>;
            return Detail::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localScvfIdx]);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Box scvf geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " type=" << type);
    }

    //! Create the sub control volume face geometries on the boundary
    ScvfCornerStorage getBoundaryScvfCorners(const Intersection& is,
                                             const typename Intersection::Geometry& isGeometry,
                                             unsigned int indexInIntersection) const
    {
        const auto refElement = referenceElement(geo_);

        const auto vIdxLocal = refElement.subEntity(is.indexInInside(), 1, indexInIntersection, dim);
        if (indexInIntersection == 0)
            return ScvfCornerStorage({ geo_.global(refElement.position(vIdxLocal, dim)), isGeometry.center() });
        else if (indexInIntersection == 1)
            return ScvfCornerStorage({ isGeometry.center(), geo_.global(refElement.position(vIdxLocal, dim)) });
        else
            DUNE_THROW(Dune::InvalidStateException, "local index exceeds the number of corners of 2d intersections");
    }

    //! get scvf normal vector for dim == 2, dimworld == 3
    template <int w = dimWorld>
    typename std::enable_if<w == 3, GlobalPosition>::type
    normal(const ScvfCornerStorage& scvfCorners,
           const std::vector<unsigned int>& scvIndices) const
    {
        const auto v1 = geo_.corner(1) - geo_.corner(0);
        const auto v2 = geo_.corner(2) - geo_.corner(0);
        const auto v3 = Dumux::crossProduct(v1, v2);
        const auto t = scvfCorners[1] - scvfCorners[0];
        GlobalPosition normal = Dumux::crossProduct(v3, t);
        normal /= normal.two_norm();

        //! ensure the right direction of the normal
        const auto v = geo_.corner(scvIndices[1]) - geo_.corner(scvIndices[0]);
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
        const auto v = geo_.corner(scvIndices[1]) - geo_.corner(scvIndices[0]);
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

    //! the wrapped element geometry
    const typename Element::Geometry& elementGeometry() const
    { return geo_; }

private:
    const typename Element::Geometry& geo_; //!< Reference to the element geometry
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

public:
    BoxGeometryHelper(const typename Element::Geometry& geometry)
    : geo_(geometry)
    {}

    //! Create a vector with the scv corners
    ScvCornerStorage getScvCorners(unsigned int localScvIdx) const
    {
        const auto type = geo_.type();
        if (type == Dune::GeometryTypes::tetrahedron)
        {
            using Corners = Detail::ScvCorners<Dune::GeometryTypes::tetrahedron>;
            return Detail::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localScvIdx]);
        }
        else if (type == Dune::GeometryTypes::prism)
        {
            using Corners = Detail::ScvCorners<Dune::GeometryTypes::prism>;
            return Detail::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localScvIdx]);
        }
        else if (type == Dune::GeometryTypes::hexahedron)
        {
            using Corners = Detail::ScvCorners<Dune::GeometryTypes::hexahedron>;
            return Detail::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localScvIdx]);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Box scv geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " type=" << type);
    }

    //! Create a vector with the scvf corners
    ScvfCornerStorage getScvfCorners(unsigned int localScvfIdx) const
    {
        // proceed according to number of corners
        const auto type = geo_.type();
        if (type == Dune::GeometryTypes::tetrahedron)
        {
            using Corners = Detail::ScvfCorners<Dune::GeometryTypes::tetrahedron>;
            return Detail::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localScvfIdx]);
        }
        else if (type == Dune::GeometryTypes::prism)
        {
            using Corners = Detail::ScvfCorners<Dune::GeometryTypes::prism>;
            return Detail::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localScvfIdx]);
        }
        else if (type == Dune::GeometryTypes::hexahedron)
        {
            using Corners = Detail::ScvfCorners<Dune::GeometryTypes::hexahedron>;
            return Detail::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localScvfIdx]);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Box scvf geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " type=" << type);
    }

    //! Create the sub control volume face geometries on the boundary
    ScvfCornerStorage getBoundaryScvfCorners(const Intersection& is,
                                             const typename Intersection::Geometry& isGeometry,
                                             unsigned int indexInIntersection) const
    {
        const auto refElement = referenceElement(geo_);
        const auto faceRefElem = referenceElement(isGeometry);

        GlobalPosition pi[9];
        auto corners = isGeometry.corners();

        // the facet center
        pi[0] = isGeometry.center();

        // corners
        const auto idxInInside = is.indexInInside();
        for (int i = 0; i < corners; ++i)
        {
            const auto vIdxLocal = refElement.subEntity(idxInInside, 1, i, dim);
            pi[i+1] = geo_.corner(vIdxLocal);
        }

        // edge midpoints
        for (int i = 0; i < faceRefElem.size(1); ++i)
        {
            const auto edgeIdxLocal = refElement.subEntity(idxInInside, 1, i, dim-1);
            pi[i+corners+1] = geo_.global(refElement.position(edgeIdxLocal, dim-1));
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

        const auto v = geo_.corner(scvIndices[1]) - geo_.corner(scvIndices[0]);
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

    //! the wrapped element geometry
    const typename Element::Geometry& elementGeometry() const
    { return geo_; }

private:
    const typename Element::Geometry& geo_; //!< Reference to the element geometry
};

} // end namespace Dumux

#endif
