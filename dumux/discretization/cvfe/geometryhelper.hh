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
 * \ingroup CvfeDiscretization
 * \brief Helper class constructing the dual grid finite volume geometries
 *        for the cvfe discretizazion method
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_CVFE_GEOMETRY_HELPER_HH

#include <array>

#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/math.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>

namespace Dumux {

template <class ct>
using CvfeMLGeometryTraits = Dune::MultiLinearGeometryTraits<ct>;


namespace Detail::Cvfe {

template<Dune::GeometryType::Id gt>
struct ScvCorners
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static const std::array<std::array<Key, 4>, 5> keys();
};

template<>
struct ScvCorners<Dune::GeometryTypes::line>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static const std::array<std::array<Key, 2>, 2> keys()
    {
        return
        {{
            { Key{0, 1}, Key{0, 0} },
            { Key{1, 1}, Key{0, 0} }
        }};
    }
};

template<>
struct ScvCorners<Dune::GeometryTypes::triangle>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static const std::array<std::vector<Key>, 4> keys()
    {
        return
        {{
            { Key{0, 2}, Key{0, 1}, Key{1, 1}, Key{0, 0} },
            { Key{1, 2}, Key{2, 1}, Key{0, 1}, Key{0, 0} },
            { Key{2, 2}, Key{1, 1}, Key{2, 1}, Key{0, 0} },
            { Key{0, 1}, Key{1, 1}, Key{2, 1} }
        }};
    }
};

template<>
struct ScvCorners<Dune::GeometryTypes::quadrilateral>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static const std::array<std::vector<Key>, 5> keys()
    {
        return
        {{
            { Key{0, 2}, Key{2, 1}, Key{0, 1}, Key{0, 0} },
            { Key{1, 2}, Key{1, 1}, Key{2, 1}, Key{0, 0} },
            { Key{2, 2}, Key{0, 1}, Key{3, 1}, Key{0, 0} },
            { Key{3, 2}, Key{3, 1}, Key{1, 1}, Key{0, 0} },
            { Key{2, 1}, Key{1, 1}, Key{0, 1}, Key{3, 1} }
        }};
    }
};

template<>
struct ScvCorners<Dune::GeometryTypes::tetrahedron>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static const std::array<std::vector<Key>, 5> keys()
    {
        return
        {{
            { Key{0, 3}, Key{0, 2}, Key{1, 2}, Key{0, 1}, Key{3, 2}, Key{1, 1}, Key{2, 1}, Key{0, 0} },
            { Key{1, 3}, Key{2, 2}, Key{0, 2}, Key{0, 1}, Key{4, 2}, Key{3, 1}, Key{1, 1}, Key{0, 0} },
            { Key{2, 3}, Key{1, 2}, Key{2, 2}, Key{0, 1}, Key{5, 2}, Key{2, 1}, Key{3, 1}, Key{0, 0} },
            { Key{3, 3}, Key{3, 2}, Key{5, 2}, Key{2, 1}, Key{4, 2}, Key{1, 1}, Key{3, 1}, Key{0, 0} },
            { Key{0, 1}, Key{1, 1}, Key{2, 1}, Key{3, 1} }
        }};
    }
};

template<>
struct ScvCorners<Dune::GeometryTypes::hexahedron>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static const std::array<std::vector<Key>, 9> keys()
    {
        return
        {{
        { Key{0, 3}, Key{6, 2}, Key{4, 2}, Key{4, 1}, Key{0, 2}, Key{2, 1}, Key{0, 1}, Key{0, 0} },
        { Key{1, 3}, Key{5, 2}, Key{6, 2}, Key{4, 1}, Key{1, 2}, Key{1, 1}, Key{2, 1}, Key{0, 0} },
        { Key{2, 3}, Key{4, 2}, Key{7, 2}, Key{4, 1}, Key{2, 2}, Key{0, 1}, Key{3, 1}, Key{0, 0} },
        { Key{3, 3}, Key{7, 2}, Key{5, 2}, Key{4, 1}, Key{3, 2}, Key{3, 1}, Key{1, 1}, Key{0, 0} },
        { Key{4, 3}, Key{8, 2}, Key{10, 2}, Key{5, 1}, Key{0, 2}, Key{0, 1}, Key{2, 1}, Key{0, 0} },
        { Key{5, 3}, Key{10, 2}, Key{9, 2}, Key{5, 1}, Key{1, 2}, Key{2, 1}, Key{1, 1}, Key{0, 0} },
        { Key{6, 3}, Key{11, 2}, Key{8, 2}, Key{5, 1}, Key{2, 2}, Key{3, 1}, Key{0, 1}, Key{0, 0} },
        { Key{7, 3}, Key{9, 2}, Key{11, 2}, Key{5, 1}, Key{3, 2}, Key{1, 1}, Key{3, 1}, Key{0, 0} },
        { Key{0, 1}, Key{2, 1}, Key{3, 1}, Key{1, 1}, Key{4, 1}, Key{5, 1} }
        }};
    }
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
    static constexpr std::array<std::array<Key, 2>, 6> keys = {{
        { Key{0, 0}, Key{0, 1} },
        { Key{1, 1}, Key{0, 0} },
        { Key{0, 0}, Key{2, 1} },
        { Key{0, 1}, Key{1, 1} },
        { Key{0, 1}, Key{2, 1} },
        { Key{1, 1}, Key{2, 1} }
    }};
};

template<>
struct ScvfCorners<Dune::GeometryTypes::quadrilateral>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 2>, 8> keys = {{
        { Key{0, 1}, Key{0, 0} },
        { Key{0, 0}, Key{1, 1} },
        { Key{0, 0}, Key{2, 1} },
        { Key{3, 1}, Key{0, 0} },
        { Key{0, 1}, Key{2, 1} },
        { Key{2, 1}, Key{1, 1} },
        { Key{0, 1}, Key{3, 1} },
        { Key{1, 1}, Key{3, 1} }
    }};
};

template<>
struct ScvfCorners<Dune::GeometryTypes::tetrahedron>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 4>, 10> keys = {{
        { Key{0, 2}, Key{0, 1}, Key{1, 1}, Key{0, 0} },
        { Key{0, 1}, Key{1, 2}, Key{0, 0}, Key{2, 1} },
        { Key{2, 2}, Key{0, 1}, Key{3, 1}, Key{0, 0} },
        { Key{2, 1}, Key{3, 2}, Key{0, 0}, Key{1, 1} },
        { Key{3, 1}, Key{0, 0}, Key{4, 2}, Key{1, 1} },
        { Key{5, 2}, Key{2, 1}, Key{3, 1}, Key{0, 0} },
        { Key{0, 1}, Key{1, 1}, Key{2, 1} },
        { Key{0, 1}, Key{1, 1}, Key{3, 1} },
        { Key{0, 1}, Key{2, 1}, Key{3, 1} },
        { Key{1, 1}, Key{2, 1}, Key{3, 1} }
    }};
};

template<>
struct ScvfCorners<Dune::GeometryTypes::hexahedron>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 4>, 20> keys = {{
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
        { Key{11, 2}, Key{5, 1}, Key{3, 1}, Key{0, 0} },
        { Key{4, 1}, Key{0, 1}, Key{2, 1} },
        { Key{4, 1}, Key{2, 1}, Key{1, 1} },
        { Key{4, 1}, Key{0, 1}, Key{3, 1} },
        { Key{4, 1}, Key{1, 1}, Key{3, 1} },
        { Key{5, 1}, Key{0, 1}, Key{2, 1} },
        { Key{5, 1}, Key{2, 1}, Key{1, 1} },
        { Key{5, 1}, Key{0, 1}, Key{3, 1} },
        { Key{5, 1}, Key{1, 1}, Key{3, 1} }
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
template<class S, class Geo, class T>
S keyToCornerStorage(const Geo& geo, const std::vector<T>& key)
{
    using Dune::referenceElement;
    const auto ref = referenceElement(geo);

    S storage(key.size());
    for(int i=0; i<key.size(); i++)
        storage[i] = geo.global(ref.position(key[i].first, key[i].second));

    return storage;
}

} // end namespace Detail

//! A class to create sub control volume and sub control volume face geometries per element
template <class GridView, class ScvType, class ScvfType>
class CvfeGeometryHelper
{
    using Scalar = typename GridView::ctype;
    using GlobalPosition = typename Dune::FieldVector<Scalar, GridView::dimensionworld>;
    using ScvCornerStorage = typename ScvType::Traits::CornerStorage;
    using ScvfCornerStorage = typename ScvfType::Traits::CornerStorage;
    using LocalIndexType = typename ScvType::Traits::LocalIndexType;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    static constexpr auto dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;
public:

    CvfeGeometryHelper(const typename Element::Geometry& geometry)
    : geo_(geometry)
    {}

    //! Create a vector with the scv corners
    ScvCornerStorage getScvCorners(unsigned int localScvIdx) const
    {
        // proceed according to number of corners of the element
        const auto type = geo_.type();
        if (type == Dune::GeometryTypes::triangle)
        {
            using Corners = Detail::Cvfe::ScvCorners<Dune::GeometryTypes::triangle>;
            return Detail::Cvfe::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys()[localScvIdx]);
        }
        else if (type == Dune::GeometryTypes::quadrilateral)
        {
            using Corners = Detail::Cvfe::ScvCorners<Dune::GeometryTypes::quadrilateral>;
            return Detail::Cvfe::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys()[localScvIdx]);
        }
        else if (type == Dune::GeometryTypes::tetrahedron)
        {
            using Corners = Detail::Cvfe::ScvCorners<Dune::GeometryTypes::tetrahedron>;
            return Detail::Cvfe::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys()[localScvIdx]);
        }
        else if (type == Dune::GeometryTypes::hexahedron)
        {
            using Corners = Detail::Cvfe::ScvCorners<Dune::GeometryTypes::hexahedron>;
            return Detail::Cvfe::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys()[localScvIdx]);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Cvfe scv geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " type=" << type);
    }

    Dune::GeometryType getScvGeometryType(unsigned int localScvIdx) const
    {
        // proceed according to number of corners of the element
        const auto type = geo_.type();
        const auto numCorners = geo_.corners();
        if (type == Dune::GeometryTypes::simplex(dim))
            return (localScvIdx < numCorners) ? Dune::GeometryTypes::cube(dim) : Dune::GeometryTypes::simplex(dim);
        else if (type == Dune::GeometryTypes::quadrilateral)
            Dune::GeometryTypes::cube(dim);
        else
            DUNE_THROW(Dune::NotImplemented, "Cvfe scv geometries for dim=" << dim
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
            using Corners = Detail::Cvfe::ScvfCorners<Dune::GeometryTypes::triangle>;
            return Detail::Cvfe::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localScvfIdx]);
        }
        else if (type == Dune::GeometryTypes::quadrilateral)
        {
            using Corners = Detail::Cvfe::ScvfCorners<Dune::GeometryTypes::quadrilateral>;
            return Detail::Cvfe::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localScvfIdx]);
        }
        else if (type == Dune::GeometryTypes::tetrahedron)
        {
            using Corners = Detail::Cvfe::ScvfCorners<Dune::GeometryTypes::tetrahedron>;
            return Detail::Cvfe::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localScvfIdx]);
        }
        else if (type == Dune::GeometryTypes::hexahedron)
        {
            using Corners = Detail::Cvfe::ScvfCorners<Dune::GeometryTypes::hexahedron>;
            return Detail::Cvfe::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localScvfIdx]);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Cvfe scvf geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " type=" << type);
    }

    Dune::GeometryType getScvfGeometryType(unsigned int localScvfIdx) const
    {
        const auto type = geo_.type();
        if (type == Dune::GeometryTypes::simplex(dim))
            return (localScvfIdx < referenceElement(geo_).size(dim-1)) ? Dune::GeometryTypes::cube(dim-1) : Dune::GeometryTypes::simplex(dim-1);
        else if (type == Dune::GeometryTypes::cube(dim))
            return (localScvfIdx < referenceElement(geo_).size(dim-1)) ? Dune::GeometryTypes::cube(dim-1) : Dune::GeometryTypes::simplex(dim-1);
        else
            DUNE_THROW(Dune::NotImplemented, "Cvfe scvf geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " type=" << type);
    }

    Dune::GeometryType getBoundaryScvfGeometryType(unsigned int localScvfIdx) const
    {
        return Dune::GeometryTypes::cube(dim-1);
    }

    //! Create the sub control volume face geometries on the boundary
    ScvfCornerStorage getBoundaryScvfCorners(unsigned int localFacetIndex,
                                             unsigned int indexInFacet) const
    {
        // We do the same as for the Box method (Remark we do use the Detail namespace from the box helper)
        if constexpr (dim == 2)
        {
            using Corners = Detail::ScvCorners<Dune::GeometryTypes::line>;
            constexpr int facetCodim = 1;
            return Detail::subEntityKeyToCornerStorage<ScvfCornerStorage>(geo_, localFacetIndex, facetCodim, Corners::keys[indexInFacet]);
        }
        else if constexpr (dim == 3)
        {
            constexpr int facetCodim = 1;
            using Dune::referenceElement;
            const auto type = referenceElement(geo_).type(localFacetIndex, facetCodim);
            if (type == Dune::GeometryTypes::triangle)
            {
                using Corners = Detail::ScvCorners<Dune::GeometryTypes::triangle>;
                return Detail::subEntityKeyToCornerStorage<ScvfCornerStorage>(geo_, localFacetIndex, facetCodim, Corners::keys[indexInFacet]);
            }
            else if (type == Dune::GeometryTypes::quadrilateral)
            {
                using Corners = Detail::ScvCorners<Dune::GeometryTypes::quadrilateral>;
                return Detail::subEntityKeyToCornerStorage<ScvfCornerStorage>(geo_, localFacetIndex, facetCodim, Corners::keys[indexInFacet]);
            }
            else
                DUNE_THROW(Dune::NotImplemented, "Cvfe boundary scvf corners for dim=" << dim
                                                                << " dimWorld=" << dimWorld
                                                                << " type=" << type);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Cvfe boundary scvf corners for dim=" << dim
                                                            << " dimWorld=" << dimWorld);
    }

    template<int d = dimWorld, std::enable_if_t<(d==3), int> = 0>
    GlobalPosition normal(const ScvfCornerStorage& p, const std::array<LocalIndexType, 2>& scvPair)
    {
        auto normal = Dumux::crossProduct(p[1]-p[0], p[2]-p[0]);
        normal /= normal.two_norm();

        GlobalPosition v = dofPosition(scvPair[1]) - dofPosition(scvPair[0]);

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

        GlobalPosition v = dofPosition(scvPair[1]) - dofPosition(scvPair[0]);

        const auto s = v*normal;
        if (std::signbit(s))
            normal *= -1;

        return normal;
    }

    //! the wrapped element geometry
    const typename Element::Geometry& elementGeometry() const
    { return geo_; }

    //! number of interior sub control volume faces
    std::size_t numInteriorScvf() const
    {
        return referenceElement(geo_).size(dim-1) + referenceElement(geo_).size(dim);
    }

    //! number of sub control volumes (number of codim-1 entities)
    std::size_t numScv() const
    {
        return referenceElement(geo_).size(dim) + 1;
    }

    //! get scv volume
    Scalar scvVolume(unsigned int localScvIdx, const ScvCornerStorage& p) const
    {
        const auto type = geo_.type();
        if (type == Dune::GeometryTypes::simplex(dim))
            return (localScvIdx < numScv()-1) ?  cubeVolume(p) : simplexVolume(p);
        else if (type == Dune::GeometryTypes::cube(dim))
            return (localScvIdx < numScv()-1) ?  cubeVolume(p) : rhombusVolume(p);
        else
            DUNE_THROW(Dune::NotImplemented, "Cvfe scv volume for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " type=" << type);
    }

    GlobalPosition scvCenter(unsigned int localScvIdx, const ScvCornerStorage& p) const
    {
        GlobalPosition center(0.0);
        // compute center point
        for (const auto& corner : p)
            center += corner;
        center /= p.size();

        return center;
    }

    //! number of boundary sub control volume faces for intersection
    template<class Is>
    std::size_t numBoundaryIsScvf(const Is& is) const
    {
        return is.geometry().corners();
    }

    template<class DofMapper>
    auto dofIndex(const DofMapper& dofMapper, const Element& element, unsigned int localScvIdx) const
    {
        if (localScvIdx < numScv()-1)
            return dofMapper.subIndex(element, localScvIdx, dim);
        else
            return dofMapper.index(element);
    }

    GlobalPosition dofPosition(unsigned int localScvIdx) const
    {
        if (localScvIdx < numScv()-1)
            return geo_.corner(localScvIdx);
        else
            return geo_.center();
    }

    std::array<LocalIndexType, 2> getScvPairForScvf(unsigned int localScvfIndex) const
    {
        const auto numEdges = referenceElement(geo_).size(dim-1);
        if(localScvfIndex < numEdges)
            return std::array{static_cast<LocalIndexType>(referenceElement(geo_).subEntity(localScvfIndex, dim-1, 0, dim)),
                              static_cast<LocalIndexType>(referenceElement(geo_).subEntity(localScvfIndex, dim-1, 1, dim))};
        else
            return std::array{static_cast<LocalIndexType>(numScv()-1),
                              static_cast<LocalIndexType>(localScvfIndex-numEdges)};
    }

    template<class Is>
    std::array<LocalIndexType, 2> getScvPairForBoundaryScvf(const Is& is, unsigned int localIsScvfIndex) const
    {
        const LocalIndexType insideScvIdx
            = static_cast<LocalIndexType>(referenceElement(geo_).subEntity(is.indexInInside(), 1, localIsScvfIndex, dim));
        return std::array{insideScvIdx, insideScvIdx};
    }

    bool isOverlapping(unsigned int localScvfIndex) const
    {
        if(localScvfIndex < referenceElement(geo_).size(dim-1))
            return false;
        else
            return true;
    }

private:
    Scalar cubeVolume(const ScvCornerStorage& p) const
    {
        if constexpr (dim == 2)
        {
            //! make sure we are using positive volumes
            //! Cross product of diagonals might be negative, depending on element orientation
            using std::abs;
            return 0.5*abs(Dumux::crossProduct(p[3]-p[0], p[2]-p[1]));
        }
        else if constexpr (dim == 3)
        {
            // after Grandy 1997, Efficient computation of volume of hexahedron
            const auto v = p[7]-p[0];
            return 1.0/6.0 * ( Dumux::tripleProduct(v, p[1]-p[0], p[3]-p[5])
                            + Dumux::tripleProduct(v, p[4]-p[0], p[5]-p[6])
                            + Dumux::tripleProduct(v, p[2]-p[0], p[6]-p[3]));
        }
    }

    Scalar simplexVolume(const ScvCornerStorage& p) const
    {
        if constexpr (dim == 2)
        {
            //! make sure we are using positive volumes
            //! Cross product of diagonals might be negative, depending on element orientation
            using std::abs;
            return 0.5*abs(Dumux::crossProduct(p[1]-p[0], p[2]-p[0]));
        }
        else if constexpr (dim == 3)
        {
            // after Grandy 1997, Efficient computation of volume of hexahedron
            const auto v = p[3]-p[0];
            using std::abs;
            return 1.0/6.0 * abs( Dumux::tripleProduct(v, p[1]-p[0], p[2]-p[0]));
        }
    }

    Scalar rhombusVolume(const ScvCornerStorage& p) const
    {
        if constexpr (dim == 2)
        {
            return cubeVolume(p);
        }
        else if constexpr (dim == 3)
        {
            // after Grandy 1997, Efficient computation of volume of hexahedron
            const auto v1 = p[4]-p[0];
            const auto v2 = p[5]-p[0];
            using std::abs;
            return 1.0/3.0 * (  abs( Dumux::tripleProduct(v1, p[1]-p[0], p[2]-p[0]))
                              + abs( Dumux::tripleProduct(v2, p[1]-p[0], p[2]-p[0])));
        }
    }

    const typename Element::Geometry& geo_; //!< Reference to the element geometry
};

} // end namespace Dumux

#endif
