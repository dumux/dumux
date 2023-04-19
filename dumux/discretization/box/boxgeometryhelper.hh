// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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

namespace Detail::Box {

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
    // key is a pair of a local sub-entity index (first) and the sub-entity's codim (second)
    return { geo.global(ref.position(key[I].first, key[I].second))... };
}

// convert key array to global corner storage
template<class S, class Geo, class T, std::size_t N, class Indices = std::make_index_sequence<N>>
S keyToCornerStorage(const Geo& geo, const std::array<T, N>& key)
{
    return keyToCornerStorageImpl<S>(geo, key, Indices{});
}

// convert key array to global corner storage
// for the i-th sub-entity of codim c (e.g. the i-th facet/codim-1-entity for boundaries)
template<class S, class Geo, class KeyArray, std::size_t... I>
S subEntityKeyToCornerStorageImpl(const Geo& geo, unsigned int i, unsigned int c, const KeyArray& key, std::index_sequence<I...>)
{
    using Dune::referenceElement;
    const auto ref = referenceElement(geo);
    // subEntity gives the subEntity number with respect to the codim-0 reference element
    // key is a pair of a local sub-entity index (first) and the sub-entity's codim (second) but here w.r.t. the sub-entity i/c
    return { geo.global(ref.position(ref.subEntity(i, c, key[I].first, c+key[I].second), c+key[I].second))... };
}

// convert key array to global corner storage
// for the i-th sub-entity of codim c (e.g. the i-th facet/codim-1-entity for boundaries)
template<class S, class Geo, class T, std::size_t N, class Indices = std::make_index_sequence<N>>
S subEntityKeyToCornerStorage(const Geo& geo, unsigned int i, unsigned int c, const std::array<T, N>& key)
{
    return subEntityKeyToCornerStorageImpl<S>(geo, i, c, key, Indices{});
}

} // end namespace Detail::Box

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

    static constexpr int dim = 1;
public:

    BoxGeometryHelper(const typename Element::Geometry& geometry)
    : geo_(geometry)
    {}

    //! Create a vector with the scv corners
    ScvCornerStorage getScvCorners(unsigned int localScvIdx) const
    {
        using Corners = Detail::Box::ScvCorners<Dune::GeometryTypes::line>;
        return Detail::Box::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localScvIdx]);
    }

    ScvGeometry scvGeometry(unsigned int localScvIdx) const
    {
        return { Dune::GeometryTypes::line, getScvCorners(localScvIdx) };
    }

    //! Create a vector with the corners of sub control volume faces
    ScvfCornerStorage getScvfCorners(unsigned int localScvfIdx) const
    {
        using Corners = Detail::Box::ScvfCorners<Dune::GeometryTypes::line>;
        return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localScvfIdx]);
    }

    //! Create the sub control volume face geometries on the boundary
    ScvfCornerStorage getBoundaryScvfCorners(unsigned int localFacetIndex,
                                             unsigned int) const
    {
        return ScvfCornerStorage{{ geo_.corner(localFacetIndex) }};
    }

    //! get scvf normal vector
    GlobalPosition normal(const ScvfCornerStorage& scvfCorners,
                          const std::vector<unsigned int>& scvIndices) const
    {
        auto normal = geo_.corner(1) - geo_.corner(0);
        normal /= normal.two_norm();
        return normal;
    }

    //! number of sub control volume faces (number of edges)
    std::size_t numInteriorScvf() const
    {
        return referenceElement(geo_).size(dim-1);
    }

    //! number of sub control volumes (number of vertices)
    std::size_t numScv() const
    {
        return referenceElement(geo_).size(dim);
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
            using Corners = Detail::Box::ScvCorners<Dune::GeometryTypes::triangle>;
            return Detail::Box::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localScvIdx]);
        }
        else if (type == Dune::GeometryTypes::quadrilateral)
        {
            using Corners = Detail::Box::ScvCorners<Dune::GeometryTypes::quadrilateral>;
            return Detail::Box::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localScvIdx]);
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
            using Corners = Detail::Box::ScvfCorners<Dune::GeometryTypes::triangle>;
            return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localScvfIdx]);
        }
        else if (type == Dune::GeometryTypes::quadrilateral)
        {
            using Corners = Detail::Box::ScvfCorners<Dune::GeometryTypes::quadrilateral>;
            return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localScvfIdx]);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Box scvf geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " type=" << type);
    }

    //! Create the sub control volume face geometries on the boundary
    ScvfCornerStorage getBoundaryScvfCorners(unsigned int localFacetIndex,
                                             unsigned int indexInFacet) const
    {
        // we have to use the corresponding facet geometry as the intersection geometry
        // might be rotated or flipped. This makes sure that the corners (dof location)
        // and corresponding scvfs are sorted in the same way
        using Corners = Detail::Box::ScvCorners<Dune::GeometryTypes::line>;
        constexpr int facetCodim = 1;
        return Detail::Box::subEntityKeyToCornerStorage<ScvfCornerStorage>(geo_, localFacetIndex, facetCodim, Corners::keys[indexInFacet]);
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

    //! number of sub control volume faces (number of edges)
    std::size_t numInteriorScvf() const
    {
        return referenceElement(geo_).size(dim-1);
    }

    //! number of sub control volumes (number of vertices)
    std::size_t numScv() const
    {
        return referenceElement(geo_).size(dim);
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
            using Corners = Detail::Box::ScvCorners<Dune::GeometryTypes::tetrahedron>;
            return Detail::Box::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localScvIdx]);
        }
        else if (type == Dune::GeometryTypes::prism)
        {
            using Corners = Detail::Box::ScvCorners<Dune::GeometryTypes::prism>;
            return Detail::Box::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localScvIdx]);
        }
        else if (type == Dune::GeometryTypes::hexahedron)
        {
            using Corners = Detail::Box::ScvCorners<Dune::GeometryTypes::hexahedron>;
            return Detail::Box::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localScvIdx]);
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
            using Corners = Detail::Box::ScvfCorners<Dune::GeometryTypes::tetrahedron>;
            return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localScvfIdx]);
        }
        else if (type == Dune::GeometryTypes::prism)
        {
            using Corners = Detail::Box::ScvfCorners<Dune::GeometryTypes::prism>;
            return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localScvfIdx]);
        }
        else if (type == Dune::GeometryTypes::hexahedron)
        {
            using Corners = Detail::Box::ScvfCorners<Dune::GeometryTypes::hexahedron>;
            return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localScvfIdx]);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Box scvf geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " type=" << type);
    }

    //! Create the sub control volume face geometries on the boundary
    ScvfCornerStorage getBoundaryScvfCorners(unsigned localFacetIndex,
                                             unsigned int indexInFacet) const
    {
        constexpr int facetCodim = 1;

        // we have to use the corresponding facet geometry as the intersection geometry
        // might be rotated or flipped. This makes sure that the corners (dof location)
        // and corresponding scvfs are sorted in the same way
        using Dune::referenceElement;
        const auto type = referenceElement(geo_).type(localFacetIndex, facetCodim);
        if (type == Dune::GeometryTypes::triangle)
        {
            using Corners = Detail::Box::ScvCorners<Dune::GeometryTypes::triangle>;
            return Detail::Box::subEntityKeyToCornerStorage<ScvfCornerStorage>(geo_, localFacetIndex, facetCodim, Corners::keys[indexInFacet]);
        }
        else if (type == Dune::GeometryTypes::quadrilateral)
        {
            using Corners = Detail::Box::ScvCorners<Dune::GeometryTypes::quadrilateral>;
            return Detail::Box::subEntityKeyToCornerStorage<ScvfCornerStorage>(geo_, localFacetIndex, facetCodim, Corners::keys[indexInFacet]);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Box boundary scvf geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " type=" << type);
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

    //! number of sub control volume faces (number of edges)
    std::size_t numInteriorScvf() const
    {
        return referenceElement(geo_).size(dim-1);
    }

    //! number of sub control volumes (number of vertices)
    std::size_t numScv() const
    {
        return referenceElement(geo_).size(dim);
    }

    //! the wrapped element geometry
    const typename Element::Geometry& elementGeometry() const
    { return geo_; }

private:
    const typename Element::Geometry& geo_; //!< Reference to the element geometry
};

} // end namespace Dumux

#endif
