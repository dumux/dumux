// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup DiamondDiscretization
 * \copydoc Dumux::DiamondGeometryHelper
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_DIAMOND_GEOMETRY_HELPER_HH

#include <array>

#include <dune/common/reservedvector.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/type.hh>

#include <dumux/common/math.hh>

namespace Dumux {

//! Traits for an efficient corner storage for fc diamond method
template <class ct>
struct FCDiamondMLGeometryTraits : public Dune::MultiLinearGeometryTraits<ct>
{
    // we use static vectors to store the corners as we know
    // the maximum number of corners in advance (2^dim)
    template< int mydim, int cdim >
    struct CornerStorage
    {
        using Type = Dune::ReservedVector< Dune::FieldVector< ct, cdim >, (1<<mydim)>;
    };
};

namespace Detail::FCDiamond {

template<Dune::GeometryType::Id gt>
struct ScvCorners;

template<>
struct ScvCorners<Dune::GeometryTypes::triangle>
{
    static constexpr Dune::GeometryType type()
    { return Dune::GeometryTypes::triangle; }

    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 3>, 3> keys = {{
        { Key{0, 2}, Key{1, 2}, Key{0, 0} },
        { Key{2, 2}, Key{0, 2}, Key{0, 0} },
        { Key{1, 2}, Key{2, 2}, Key{0, 0} }
    }};
};

template<>
struct ScvCorners<Dune::GeometryTypes::quadrilateral>
{
    static constexpr Dune::GeometryType type()
    { return Dune::GeometryTypes::triangle; }

    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 3>, 4> keys = {{
        { Key{2, 2}, Key{0, 2}, Key{0, 0} },
        { Key{1, 2}, Key{3, 2}, Key{0, 0} },
        { Key{0, 2}, Key{1, 2}, Key{0, 0} },
        { Key{3, 2}, Key{2, 2}, Key{0, 0} }
    }};
};

template<>
struct ScvCorners<Dune::GeometryTypes::tetrahedron>
{
    static constexpr Dune::GeometryType type()
    { return Dune::GeometryTypes::tetrahedron; }

    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 4>, 4> keys = {{
        { Key{0, 3}, Key{1, 3}, Key{2, 3}, Key{0, 0} },
        { Key{1, 3}, Key{0, 3}, Key{3, 3}, Key{0, 0} },
        { Key{0, 3}, Key{2, 3}, Key{3, 3}, Key{0, 0} },
        { Key{2, 3}, Key{1, 3}, Key{3, 3}, Key{0, 0} }
    }};
};

template<>
struct ScvCorners<Dune::GeometryTypes::hexahedron>
{
    static constexpr Dune::GeometryType type()
    { return Dune::GeometryTypes::pyramid; }

    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 5>, 6> keys = {{
        { Key{4, 3}, Key{0, 3}, Key{6, 3}, Key{2, 3}, Key{0, 0} },
        { Key{1, 3}, Key{5, 3}, Key{3, 3}, Key{7, 3}, Key{0, 0} },
        { Key{4, 3}, Key{5, 3}, Key{0, 3}, Key{1, 3}, Key{0, 0} },
        { Key{2, 3}, Key{3, 3}, Key{6, 3}, Key{7, 3}, Key{0, 0} },
        { Key{0, 3}, Key{1, 3}, Key{2, 3}, Key{3, 3}, Key{0, 0} },
        { Key{6, 3}, Key{7, 3}, Key{4, 3}, Key{5, 3}, Key{0, 0} }
    }};
};

template<Dune::GeometryType::Id gt>
struct ScvfCorners;

template<>
struct ScvfCorners<Dune::GeometryTypes::triangle>
{
    static constexpr Dune::GeometryType type()
    { return Dune::GeometryTypes::line; }

    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 2>, 3> keys = {{
        { Key{0, 2}, Key{0, 0} },
        { Key{1, 2}, Key{0, 0} },
        { Key{2, 2}, Key{0, 0} }
    }};
};

template<>
struct ScvfCorners<Dune::GeometryTypes::quadrilateral>
{
    static constexpr Dune::GeometryType type()
    { return Dune::GeometryTypes::line; }

    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 2>, 4> keys = {{
        { Key{0, 2}, Key{0, 0} },
        { Key{1, 2}, Key{0, 0} },
        { Key{2, 2}, Key{0, 0} },
        { Key{3, 2}, Key{0, 0} }
    }};
};

template<>
struct ScvfCorners<Dune::GeometryTypes::tetrahedron>
{
    static constexpr Dune::GeometryType type()
    { return Dune::GeometryTypes::triangle; }

    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 3>, 6> keys = {{
        { Key{0, 3}, Key{1, 3}, Key{0, 0} },
        { Key{2, 3}, Key{0, 3}, Key{0, 0} },
        { Key{1, 3}, Key{2, 3}, Key{0, 0} },
        { Key{0, 3}, Key{3, 3}, Key{0, 0} },
        { Key{1, 3}, Key{3, 3}, Key{0, 0} },
        { Key{2, 3}, Key{3, 3}, Key{0, 0} }
    }};
};

template<>
struct ScvfCorners<Dune::GeometryTypes::hexahedron>
{
    static constexpr Dune::GeometryType type()
    { return Dune::GeometryTypes::triangle; }

    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 3>, 12> keys = {{
        { Key{0, 3}, Key{4, 3}, Key{0, 0} },
        { Key{1, 3}, Key{5, 3}, Key{0, 0} },
        { Key{2, 3}, Key{6, 3}, Key{0, 0} },
        { Key{3, 3}, Key{7, 3}, Key{0, 0} },
        { Key{0, 3}, Key{2, 3}, Key{0, 0} },
        { Key{1, 3}, Key{3, 3}, Key{0, 0} },
        { Key{0, 3}, Key{1, 3}, Key{0, 0} },
        { Key{2, 3}, Key{3, 3}, Key{0, 0} },
        { Key{4, 3}, Key{6, 3}, Key{0, 0} },
        { Key{5, 3}, Key{7, 3}, Key{0, 0} },
        { Key{4, 3}, Key{5, 3}, Key{0, 0} },
        { Key{6, 3}, Key{7, 3}, Key{0, 0} }
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

// boundary corners for the i-th facet
template<class S, class Geo, std::size_t... ii>
S boundaryCornerStorageImpl(const Geo& geo, unsigned int i, std::index_sequence<ii...>)
{
    using Dune::referenceElement;
    const auto ref = referenceElement(geo);
    // simply the vertices of the facet i
    return { geo.global(ref.position(ref.subEntity(i, 1, ii, Geo::mydimension), Geo::mydimension))... };
}

// boundary corners for the i-th facet
template<class S, std::size_t numCorners, class Geo>
S boundaryCornerStorage(const Geo& geo, unsigned int i)
{
    return boundaryCornerStorageImpl<S>(geo, i, std::make_index_sequence<numCorners>{});
}

template<class IndexType, Dune::GeometryType::Id gt>
struct InsideOutsideScv;

template<class IndexType>
struct InsideOutsideScv<IndexType, Dune::GeometryTypes::triangle>
{
    static constexpr std::array<std::array<IndexType, 2>, 3> pairs = {{
        {0, 1}, {0, 2}, {1, 2}
    }};
};

template<class IndexType>
struct InsideOutsideScv<IndexType, Dune::GeometryTypes::quadrilateral>
{
    static constexpr std::array<std::array<IndexType, 2>, 4> pairs = {{
        {0, 2}, {1, 2}, {0, 3}, {1, 3}
    }};
};

template<class IndexType>
struct InsideOutsideScv<IndexType, Dune::GeometryTypes::tetrahedron>
{
    static constexpr std::array<std::array<IndexType, 2>, 6> pairs = {{
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3},
    }};
};

template<class IndexType>
struct InsideOutsideScv<IndexType, Dune::GeometryTypes::hexahedron>
{
    static constexpr std::array<std::array<IndexType, 2>, 12> pairs = {{
        {0, 2}, {1, 2}, {0, 3}, {1, 3}, {0, 4}, {1, 4},
        {2, 4}, {3, 4}, {0, 5}, {1, 5}, {2, 5}, {3, 5}
    }};
};

} // end namespace Detail::FCDiamond

/*!
 * \ingroup DiamondDiscretization
 * \brief Helper class to construct SCVs and SCVFs for the diamond scheme
 */
template <class GridView, class ScvType, class ScvfType>
class DiamondGeometryHelper
{
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr auto dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;

    using ScvCornerStorage = typename ScvType::Traits::CornerStorage;
    using LocalIndexType = typename ScvType::Traits::LocalIndexType;
    using ScvfCornerStorage = typename ScvfType::Traits::CornerStorage;

    // for the normal
    using ctype = typename GridView::ctype;
    using GlobalPosition = typename Dune::FieldVector<ctype, GridView::dimensionworld>;

public:
    DiamondGeometryHelper(const typename Element::Geometry& geo)
    : geo_(geo)
    {}

    //! Create a corner storage with the scv corners for a given face (codim-1) index
    ScvCornerStorage getScvCorners(unsigned int localFacetIndex) const
    {
        const auto type = geo_.type();
        if (type == Dune::GeometryTypes::triangle)
        {
            using Corners = Detail::FCDiamond::ScvCorners<Dune::GeometryTypes::triangle>;
            return Detail::FCDiamond::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localFacetIndex]);
        }
        else if (type == Dune::GeometryTypes::quadrilateral)
        {
            using Corners = Detail::FCDiamond::ScvCorners<Dune::GeometryTypes::quadrilateral>;
            return Detail::FCDiamond::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localFacetIndex]);
        }
        else if (type == Dune::GeometryTypes::tetrahedron)
        {
            using Corners = Detail::FCDiamond::ScvCorners<Dune::GeometryTypes::tetrahedron>;
            return Detail::FCDiamond::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localFacetIndex]);
        }
        else if (type == Dune::GeometryTypes::hexahedron)
        {
            using Corners = Detail::FCDiamond::ScvCorners<Dune::GeometryTypes::hexahedron>;
            return Detail::FCDiamond::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localFacetIndex]);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Scv geometries for type " << type);
    }

    //! Create a corner storage with the scvf corners for a given edge (codim-2) index
    ScvfCornerStorage getScvfCorners(unsigned int localEdgeIndex) const
    {
        const auto type = geo_.type();
        if (type == Dune::GeometryTypes::triangle)
        {
            using Corners = Detail::FCDiamond::ScvfCorners<Dune::GeometryTypes::triangle>;
            return Detail::FCDiamond::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localEdgeIndex]);
        }
        else if (type == Dune::GeometryTypes::quadrilateral)
        {
            using Corners = Detail::FCDiamond::ScvfCorners<Dune::GeometryTypes::quadrilateral>;
            return Detail::FCDiamond::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localEdgeIndex]);
        }
        else if (type == Dune::GeometryTypes::tetrahedron)
        {
            using Corners = Detail::FCDiamond::ScvfCorners<Dune::GeometryTypes::tetrahedron>;
            return Detail::FCDiamond::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localEdgeIndex]);
        }
        else if (type == Dune::GeometryTypes::hexahedron)
        {
            using Corners = Detail::FCDiamond::ScvfCorners<Dune::GeometryTypes::hexahedron>;
            return Detail::FCDiamond::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localEdgeIndex]);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Scvf geometries for type " << type);
    }

    //! Create the sub control volume face geometries on the boundary for a given face index
    ScvfCornerStorage getBoundaryScvfCorners(unsigned int localFacetIndex) const
    {
        const auto type = geo_.type();
        if (type == Dune::GeometryTypes::triangle || type == Dune::GeometryTypes::quadrilateral)
            return Detail::FCDiamond::boundaryCornerStorage<ScvfCornerStorage, 2>(geo_, localFacetIndex);
        else if (type == Dune::GeometryTypes::tetrahedron)
            return Detail::FCDiamond::boundaryCornerStorage<ScvfCornerStorage, 3>(geo_, localFacetIndex);
        else if (type == Dune::GeometryTypes::hexahedron)
            return Detail::FCDiamond::boundaryCornerStorage<ScvfCornerStorage, 4>(geo_, localFacetIndex);
        else
            DUNE_THROW(Dune::NotImplemented, "Boundary scvf geometries for type " << type);
    }

    GlobalPosition facetCenter(unsigned int localFacetIndex) const
    {
        return geo_.global(referenceElement(geo_).position(localFacetIndex, 1));
    }

    std::array<LocalIndexType, 2> getInsideOutsideScvForScvf(unsigned int localEdgeIndex)
    {
        const auto type = geo_.type();
        if (type == Dune::GeometryTypes::triangle)
            return Detail::FCDiamond::InsideOutsideScv<LocalIndexType, Dune::GeometryTypes::triangle>::pairs[localEdgeIndex];
        else if (type == Dune::GeometryTypes::quadrilateral)
            return Detail::FCDiamond::InsideOutsideScv<LocalIndexType, Dune::GeometryTypes::quadrilateral>::pairs[localEdgeIndex];
        else if (type == Dune::GeometryTypes::tetrahedron)
            return Detail::FCDiamond::InsideOutsideScv<LocalIndexType, Dune::GeometryTypes::tetrahedron>::pairs[localEdgeIndex];
        else if (type == Dune::GeometryTypes::hexahedron)
            return Detail::FCDiamond::InsideOutsideScv<LocalIndexType, Dune::GeometryTypes::hexahedron>::pairs[localEdgeIndex];
        else
            DUNE_THROW(Dune::NotImplemented, "Inside outside scv pairs for type " << type);
    }

    //! number of interior sub control volume faces (number of codim-2 entities)
    std::size_t numInteriorScvf()
    {
        return referenceElement(geo_).size(2);
    }

    //! number of sub control volumes (number of codim-1 entities)
    std::size_t numScv()
    {
        return referenceElement(geo_).size(1);
    }

    template<int d = dimWorld, std::enable_if_t<(d==3), int> = 0>
    GlobalPosition normal(const ScvfCornerStorage& p, const std::array<LocalIndexType, 2>& scvPair)
    {
        auto normal = Dumux::crossProduct(p[1]-p[0], p[2]-p[0]);
        normal /= normal.two_norm();

        const auto ref = referenceElement(geo_);
        const auto v = facetCenter_(scvPair[1], ref) - facetCenter_(scvPair[0], ref);

        const auto s = v*normal;
        if (std::signbit(s))
            normal *= -1;

        return normal;
    }

    template<int d = dimWorld, std::enable_if_t<(d==2), int> = 0>
    GlobalPosition normal(const ScvfCornerStorage& p, const std::array<LocalIndexType, 2>& scvPair)
    {
        //! obtain normal vector by 90° counter-clockwise rotation of t
        const auto t = p[1] - p[0];
        GlobalPosition normal({-t[1], t[0]});
        normal /= normal.two_norm();

        const auto ref = referenceElement(geo_);
        const auto v = facetCenter_(scvPair[1], ref) - facetCenter_(scvPair[0], ref);

        const auto s = v*normal;
        if (std::signbit(s))
            normal *= -1;

        return normal;
    }

    const typename Element::Geometry& elementGeometry() const
    { return geo_; }

private:
    template<class RefElement>
    GlobalPosition facetCenter_(unsigned int localFacetIndex, const RefElement& ref) const
    {
        return geo_.global(ref.position(localFacetIndex, 1));
    }

    const typename Element::Geometry& geo_; //!< Reference to the element geometry
};

} // end namespace Dumux

#endif
