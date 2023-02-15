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
 * \ingroup PQ2Discretization
 * \brief Helper class constructing the dual grid finite volume geometries
 *        for the cvfe discretizazion method
 */
#ifndef DUMUX_DISCRETIZATION_PQ2_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_PQ2_GEOMETRY_HELPER_HH

#include <array>

#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/common/reservedvector.hh>

#include <dumux/common/math.hh>
#include <dumux/geometry/volume.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>

namespace Dumux {

//! Traits for an efficient corner storage for the PQ2 method
template <class ct>
struct PQ2MLGeometryTraits : public Dune::MultiLinearGeometryTraits<ct>
{
    // we use static vectors to store the corners as we know
    // the maximum number of corners in advance
    template< int mydim, int cdim >
    struct CornerStorage
    {
        using Type = Dune::ReservedVector< Dune::FieldVector< ct, cdim >, (1<<mydim)+ mydim*(1<<(mydim-1))>;
    };
};

namespace Detail::PQ2 {

template<Dune::GeometryType::Id gt>
struct EdgeScvCorners;

template<>
struct EdgeScvCorners<Dune::GeometryTypes::line>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 2>, 1> keys = {{
        { Key{0, 1}, Key{1, 1} }
    }};
};

template<>
struct EdgeScvCorners<Dune::GeometryTypes::triangle>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 3>, 3> keys = {{
        { Key{0, 2}, Key{1, 2}, Key{0, 0} },
        { Key{2, 2}, Key{0, 2}, Key{0, 0} },
        { Key{1, 2}, Key{2, 2}, Key{0, 0} }
    }};
};

// template<>
// struct EdgeScvCorners<Dune::GeometryTypes::quadrilateral>
// {
//     using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
//     static constexpr std::array<std::array<Key, 4>, 1> keys = {{
//         { Key{2, 1}, Key{1, 1}, Key{0, 1}, Key{3, 1} }
//     }};
// };

// template<>
// struct EdgeScvCorners<Dune::GeometryTypes::tetrahedron>
// {
//     using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
//     static constexpr std::array<std::array<Key, 4>, 1> keys = {{
//         { Key{0, 1}, Key{1, 1}, Key{2, 1}, Key{3, 1} }
//     }};
// };

// template<>
// struct EdgeScvCorners<Dune::GeometryTypes::hexahedron>
// {
//     using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
//     static constexpr std::array<std::array<Key, 6>, 1> keys = {{
//         { Key{0, 1}, Key{2, 1}, Key{3, 1}, Key{1, 1}, Key{4, 1}, Key{5, 1} }
//     }};
// };


template<Dune::GeometryType::Id gt>
struct EdgeScvfCorners;

template<>
struct EdgeScvfCorners<Dune::GeometryTypes::triangle>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 2>, 6> keys = {{
        { Key{0, 2}, Key{0, 0} },
        { Key{1, 2}, Key{0, 0} },
        { Key{2, 2}, Key{0, 0} },
        { Key{0, 2}, Key{0, 0} },
        { Key{1, 2}, Key{0, 0} },
        { Key{2, 2}, Key{0, 0} }
    }};
};

// template<>
// struct EdgeScvfCorners<Dune::GeometryTypes::quadrilateral>
// {
//     using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
//     static constexpr std::array<std::array<Key, 2>, 4> keys = {{
//         { Key{0, 1}, Key{2, 1} },
//         { Key{2, 1}, Key{1, 1} },
//         { Key{0, 1}, Key{3, 1} },
//         { Key{1, 1}, Key{3, 1} }
//     }};
// };

// template<>
// struct EdgeScvfCorners<Dune::GeometryTypes::tetrahedron>
// {
//     using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
//     static constexpr std::array<std::array<Key, 3>, 10> keys = {{
//         { Key{0, 1}, Key{1, 1}, Key{2, 1} },
//         { Key{0, 1}, Key{1, 1}, Key{3, 1} },
//         { Key{0, 1}, Key{2, 1}, Key{3, 1} },
//         { Key{1, 1}, Key{2, 1}, Key{3, 1} }
//     }};
// };

// template<>
// struct EdgeScvfCorners<Dune::GeometryTypes::hexahedron>
// {
//     using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
//     static constexpr std::array<std::array<Key, 3>, 8> keys = {{
//         { Key{4, 1}, Key{0, 1}, Key{2, 1} },
//         { Key{4, 1}, Key{2, 1}, Key{1, 1} },
//         { Key{4, 1}, Key{0, 1}, Key{3, 1} },
//         { Key{4, 1}, Key{1, 1}, Key{3, 1} },
//         { Key{5, 1}, Key{0, 1}, Key{2, 1} },
//         { Key{5, 1}, Key{2, 1}, Key{1, 1} },
//         { Key{5, 1}, Key{0, 1}, Key{3, 1} },
//         { Key{5, 1}, Key{1, 1}, Key{3, 1} }
//     }};
// };


template<class IndexType, Dune::GeometryType::Id gt>
struct InsideOutsideOverlappingScv;

template<class IndexType>
struct InsideOutsideOverlappingScv<IndexType, Dune::GeometryTypes::triangle>
{
    static constexpr std::array<std::array<IndexType, 2>, 6> pairs = {{
        {3, 4}, {3, 5}, {4, 5}, {4, 3}, {5, 3}, {5, 4}
    }};
};


} // end namespace Detail::PQ2

/*!
 * \ingroup PQ2Discretization
 * \brief A class to create sub control volume and sub control volume face geometries per element
 */
template <class GridView, class ScvType, class ScvfType>
class PQ2GeometryHelper
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

    PQ2GeometryHelper(const typename Element::Geometry& geometry)
    : geo_(geometry)
    , boxHelper_(geometry)
    {}

    //! Create a vector with the scv corners
    ScvCornerStorage getScvCorners(unsigned int localScvIdx) const
    {
        // proceed according to number of corners of the element
        const auto type = geo_.type();
        const auto numBoxScv = boxHelper_.numScv();
        // reuse box geometry helper for the corner scvs
        if (localScvIdx < numBoxScv)
            return boxHelper_.getScvCorners(localScvIdx);

        const auto localOverlappingScvIdx = localScvIdx-numBoxScv;
        if (type == Dune::GeometryTypes::triangle)
        {
            using Corners = Detail::PQ2::EdgeScvCorners<Dune::GeometryTypes::triangle>;
            return Detail::Box::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localOverlappingScvIdx]);
        }
        // else if (type == Dune::GeometryTypes::quadrilateral)
        // {
        //     using Corners = Detail::PQ2::EdgeScvCorners<Dune::GeometryTypes::quadrilateral>;
        //     return Detail::Box::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localOverlappingScvIdx]);
        // }
        // else if (type == Dune::GeometryTypes::tetrahedron)
        // {
        //     using Corners = Detail::PQ2::EdgeScvCorners<Dune::GeometryTypes::tetrahedron>;
        //     return Detail::Box::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localOverlappingScvIdx]);
        // }
        // else if (type == Dune::GeometryTypes::hexahedron)
        // {
        //     using Corners = Detail::PQ2::EdgeScvCorners<Dune::GeometryTypes::hexahedron>;
        //     return Detail::Box::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localOverlappingScvIdx]);
        // }
        else
            DUNE_THROW(Dune::NotImplemented, "PQ2 scv geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " type=" << type);
    }

    Dune::GeometryType getScvGeometryType(unsigned int localScvIdx) const
    {
        // proceed according to number of corners of the element
        const auto type = geo_.type();
        const auto numBoxScv = boxHelper_.numScv();
        if (localScvIdx < numBoxScv)
            return Dune::GeometryTypes::cube(dim);
        else if (type == Dune::GeometryTypes::simplex(dim))
            return Dune::GeometryTypes::simplex(dim);
        else if (type == Dune::GeometryTypes::cube(dim))
            return Dune::GeometryTypes::simplex(dim);
        else
            DUNE_THROW(Dune::NotImplemented, "PQ2 scv geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " type=" << type);
    }

    //! Create a vector with the corners of sub control volume faces
    ScvfCornerStorage getScvfCorners(unsigned int localScvfIdx) const
    {
        // proceed according to number of corners
        const auto type = geo_.type();
        const auto numBoxScvf = boxHelper_.numInteriorScvf();
        // reuse box geometry helper for the corner scvs
        if (localScvfIdx < numBoxScvf)
            return boxHelper_.getScvfCorners(localScvfIdx);

        const auto localOverlappingScvfIdx = localScvfIdx-numBoxScvf;
        if (type == Dune::GeometryTypes::triangle)
        {
            using Corners = Detail::PQ2::EdgeScvfCorners<Dune::GeometryTypes::triangle>;
            return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localOverlappingScvfIdx]);
        }
        // else if (type == Dune::GeometryTypes::quadrilateral)
        // {
        //     using Corners = Detail::PQ2::EdgeScvfCorners<Dune::GeometryTypes::quadrilateral>;
        //     return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localOverlappingScvfIdx]);
        // }
        // else if (type == Dune::GeometryTypes::tetrahedron)
        // {
        //     using Corners = Detail::PQ2::EdgeScvfCorners<Dune::GeometryTypes::tetrahedron>;
        //     return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localOverlappingScvfIdx]);
        // }
        // else if (type == Dune::GeometryTypes::hexahedron)
        // {
        //     using Corners = Detail::PQ2::EdgeScvfCorners<Dune::GeometryTypes::hexahedron>;
        //     return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localOverlappingScvfIdx]);
        // }
        else
            DUNE_THROW(Dune::NotImplemented, "PQ2 scvf geometries for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " type=" << type);
    }

    Dune::GeometryType getInteriorScvfGeometryType(unsigned int localScvfIdx) const
    {
        const auto numBoxScvf = boxHelper_.numInteriorScvf();
        if (localScvfIdx < numBoxScvf)
            return Dune::GeometryTypes::cube(dim-1);
        else
            return Dune::GeometryTypes::simplex(dim-1);
    }

    //! Create the sub control volume face geometries on the boundary
    ScvfCornerStorage getBoundaryScvfCorners(unsigned int localFacetIndex,
                                             unsigned int indexInFacet) const
    {
        const auto numBoxScvf = referenceElement(geo_).size(localFacetIndex, 1, dim);
        if (indexInFacet < numBoxScvf)
        {
            return boxHelper_.getBoundaryScvfCorners(localFacetIndex, indexInFacet);
        }
        else
        {
            indexInFacet -= numBoxScvf;
            using Corners = Detail::PQ2::EdgeScvCorners<Dune::GeometryTypes::line>;
            constexpr int facetCodim = 1;
            return Detail::Box::subEntityKeyToCornerStorage<ScvfCornerStorage>(geo_, localFacetIndex, facetCodim, Corners::keys[indexInFacet]);
        }
    }

    Dune::GeometryType getBoundaryScvfGeometryType(unsigned int localFacetIndex,
                                                   unsigned int indexInFacet) const
    {
        const auto numBoxScvf = referenceElement(geo_).size(localFacetIndex, 1, dim);
        if (indexInFacet < numBoxScvf)
            return Dune::GeometryTypes::cube(dim-1);
        else
            return Dune::GeometryTypes::simplex(dim-1);
    }


    template<int d = dimWorld, std::enable_if_t<(d==3), int> = 0>
    GlobalPosition normal(const ScvfCornerStorage& p, const std::array<LocalIndexType, 2>& scvPair)
    {
        auto normal = Dumux::crossProduct(p[1]-p[0], p[2]-p[0]);
        normal /= normal.two_norm();

        auto center = std::accumulate(p.begin(), p.end(), GlobalPosition(0.0));
        center /= p.size();
        GlobalPosition v = center - dofPosition(scvPair[0]);

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

        auto center = std::accumulate(p.begin(), p.end(), GlobalPosition(0.0));
        center /= p.size();
        GlobalPosition v = center - dofPosition(scvPair[0]);

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
        return boxHelper_.numInteriorScvf() + 2*referenceElement(geo_).size(dim-1);
    }

    //! number of boundary sub control volume faces for face localFacetIndex
    std::size_t numBoundaryScvf(unsigned int localFacetIndex) const
    {
        return referenceElement(geo_).size(localFacetIndex, 1, dim) + referenceElement(geo_).size(localFacetIndex, 1, dim-1);
    }

    //! number of sub control volumes (number of codim-1 entities)
    std::size_t numScv() const
    {
        return boxHelper_.numScv() + referenceElement(geo_).size(dim-1);
    }

    //! get scv volume
    Scalar scvVolume(unsigned int localScvIdx, const ScvCornerStorage& p) const
    {
        const auto scvType = getScvGeometryType(localScvIdx);
        return Dumux::convexPolytopeVolume<dim>(
            scvType,
            [&](unsigned int i){ return p[i]; }
        );
    }

    template<class DofMapper>
    auto dofIndex(const DofMapper& dofMapper, const Element& element, unsigned int localScvIdx) const
    {
        const auto numBoxScv = boxHelper_.numScv();
        if (localScvIdx < numBoxScv)
            return dofMapper.subIndex(element, localScvIdx, dim);
        else
        {
            const auto localEdgeScvIdx = localScvIdx-numBoxScv;
            return dofMapper.subIndex(element, localEdgeScvIdx, dim-1);
        }
    }

    GlobalPosition dofPosition(unsigned int localScvIdx) const
    {
        const auto numBoxScv = boxHelper_.numScv();
        if (localScvIdx < numBoxScv)
            return geo_.corner(localScvIdx);
        else
        {
            const auto ref = referenceElement(geo_);
            const auto localEdgeScvIdx = localScvIdx-numBoxScv;
            return geo_.global(ref.position(localEdgeScvIdx, dim-1));
        }
    }

    std::array<LocalIndexType, 2> getScvPairForScvf(unsigned int localScvfIndex) const
    {
        const auto numBoxFaces = boxHelper_.numInteriorScvf();
        if (localScvfIndex < numBoxFaces)
        {
            return {
                static_cast<LocalIndexType>(referenceElement(geo_).subEntity(localScvfIndex, dim-1, 0, dim)),
                static_cast<LocalIndexType>(referenceElement(geo_).subEntity(localScvfIndex, dim-1, 1, dim))
            };
        }
        else
        {
            const auto type = geo_.type();
            if (type == Dune::GeometryTypes::triangle)
                return Detail::PQ2::InsideOutsideOverlappingScv<LocalIndexType, Dune::GeometryTypes::triangle>::pairs[localScvfIndex-numBoxFaces];

            DUNE_THROW(Dune::NotImplemented, "PQ2 getScvPairForScvf");
        }
    }

    std::array<LocalIndexType, 2> getScvPairForBoundaryScvf(unsigned int localFacetIndex, unsigned int localIsScvfIndex) const
    {
        const auto numBoxScvf = referenceElement(geo_).size(localFacetIndex, 1, dim);
        if (localIsScvfIndex < numBoxScvf)
        {
            const LocalIndexType insideScvIdx
                = static_cast<LocalIndexType>(referenceElement(geo_).subEntity(localFacetIndex, 1, localIsScvfIndex, dim));
            return { insideScvIdx, insideScvIdx };
        }
        else
        {
            localIsScvfIndex -= numBoxScvf;
            const LocalIndexType insideScvIdx
                = static_cast<LocalIndexType>(referenceElement(geo_).subEntity(localFacetIndex, 1, localIsScvfIndex, dim-1))
                + boxHelper_.numScv();
            return { insideScvIdx, insideScvIdx };
        }
    }

    bool isOverlappingBoundaryScvf(unsigned int localFacetIndex, unsigned int localIsScvfIndex) const
    {
        const auto numBoxScvf = referenceElement(geo_).size(localFacetIndex, 1, dim);
        if (localIsScvfIndex < numBoxScvf)
            return false;
        else
            return true;
    }

    bool isOverlappingScvf(unsigned int localScvfIndex) const
    {
        if (localScvfIndex < boxHelper_.numInteriorScvf())
            return false;
        else
            return true;
    }

    bool isOverlappingScv(unsigned int localScvIndex) const
    {
        if (localScvIndex < boxHelper_.numScv())
            return false;
        else
            return true;
    }

    template<class LocalKey>
    LocalIndexType localKeyToLocalScvIndex(const LocalKey& localKey) const
    {
        if(localKey.codim() == dim)
            return localKey.subEntity();
        else
            return boxHelper_.numScv() + localKey.subEntity();
    }

private:
    const typename Element::Geometry& geo_; //!< Reference to the element geometry
    Dumux::BoxGeometryHelper<GridView, dim, ScvType, ScvfType> boxHelper_;
};

} // end namespace Dumux

#endif
