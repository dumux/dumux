// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ1BubbleDiscretization
 * \brief Helper class constructing the dual grid finite volume geometries
 *        for the cvfe discretizazion method
 */
#ifndef DUMUX_DISCRETIZATION_PQ1BUBBLE_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_PQ1BUBBLE_GEOMETRY_HELPER_HH

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

//! Traits for an efficient corner storage for the PQ1Bubble method
template <class ct>
struct PQ1BubbleMLGeometryTraits : public Dune::MultiLinearGeometryTraits<ct>
{
    // we use static vectors to store the corners as we know
    // the maximum number of corners in advance (2^dim)
    template< int mydim, int cdim >
    struct CornerStorage
    {
        using Type = Dune::ReservedVector< Dune::FieldVector< ct, cdim >, (1<<mydim)+1>;
    };
};

namespace Detail::PQ1Bubble {

template<Dune::GeometryType::Id gt>
struct OverlappingScvCorners;

template<>
struct OverlappingScvCorners<Dune::GeometryTypes::line>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 2>, 1> keys = {{
        { Key{0, 1}, Key{1, 1} }
    }};
};

template<>
struct OverlappingScvCorners<Dune::GeometryTypes::triangle>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 3>, 1> keys = {{
        { Key{0, 1}, Key{1, 1}, Key{2, 1} }
    }};
};

template<>
struct OverlappingScvCorners<Dune::GeometryTypes::quadrilateral>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 4>, 1> keys = {{
        { Key{2, 1}, Key{1, 1}, Key{0, 1}, Key{3, 1} }
    }};
};

template<>
struct OverlappingScvCorners<Dune::GeometryTypes::tetrahedron>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 4>, 1> keys = {{
        { Key{0, 1}, Key{1, 1}, Key{2, 1}, Key{3, 1} }
    }};
};

template<>
struct OverlappingScvCorners<Dune::GeometryTypes::hexahedron>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 6>, 1> keys = {{
        { Key{0, 1}, Key{2, 1}, Key{3, 1}, Key{1, 1}, Key{4, 1}, Key{5, 1} }
    }};
};


template<Dune::GeometryType::Id gt>
struct OverlappingScvfCorners;

template<>
struct OverlappingScvfCorners<Dune::GeometryTypes::line>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 1>, 1> keys = {{
        { Key{0, 0} }
    }};
};

template<>
struct OverlappingScvfCorners<Dune::GeometryTypes::triangle>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 2>, 3> keys = {{
        { Key{0, 1}, Key{1, 1} },
        { Key{0, 1}, Key{2, 1} },
        { Key{1, 1}, Key{2, 1} }
    }};
};

template<>
struct OverlappingScvfCorners<Dune::GeometryTypes::quadrilateral>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 2>, 4> keys = {{
        { Key{0, 1}, Key{2, 1} },
        { Key{2, 1}, Key{1, 1} },
        { Key{0, 1}, Key{3, 1} },
        { Key{1, 1}, Key{3, 1} }
    }};
};

template<>
struct OverlappingScvfCorners<Dune::GeometryTypes::tetrahedron>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 3>, 10> keys = {{
        { Key{0, 1}, Key{1, 1}, Key{2, 1} },
        { Key{0, 1}, Key{1, 1}, Key{3, 1} },
        { Key{0, 1}, Key{2, 1}, Key{3, 1} },
        { Key{1, 1}, Key{2, 1}, Key{3, 1} }
    }};
};

template<>
struct OverlappingScvfCorners<Dune::GeometryTypes::hexahedron>
{
    using Key = std::pair<std::uint8_t, std::uint8_t>; // (i, codim)
    static constexpr std::array<std::array<Key, 3>, 8> keys = {{
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

} // end namespace Detail::PQ1Bubble

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief A class to create sub control volume and sub control volume face geometries per element
 */
template <class GridView, class ScvType, class ScvfType>
class PQ1BubbleGeometryHelper
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

    PQ1BubbleGeometryHelper(const typename Element::Geometry& geometry)
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
            using Corners = Detail::PQ1Bubble::OverlappingScvCorners<Dune::GeometryTypes::triangle>;
            return Detail::Box::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localOverlappingScvIdx]);
        }
        else if (type == Dune::GeometryTypes::quadrilateral)
        {
            using Corners = Detail::PQ1Bubble::OverlappingScvCorners<Dune::GeometryTypes::quadrilateral>;
            return Detail::Box::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localOverlappingScvIdx]);
        }
        else if (type == Dune::GeometryTypes::tetrahedron)
        {
            using Corners = Detail::PQ1Bubble::OverlappingScvCorners<Dune::GeometryTypes::tetrahedron>;
            return Detail::Box::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localOverlappingScvIdx]);
        }
        else if (type == Dune::GeometryTypes::hexahedron)
        {
            using Corners = Detail::PQ1Bubble::OverlappingScvCorners<Dune::GeometryTypes::hexahedron>;
            return Detail::Box::keyToCornerStorage<ScvCornerStorage>(geo_, Corners::keys[localOverlappingScvIdx]);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "PQ1Bubble scv geometries for dim=" << dim
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
        else if (type == Dune::GeometryTypes::quadrilateral)
            return Dune::GeometryTypes::quadrilateral;
        else if (type == Dune::GeometryTypes::hexahedron)
            return Dune::GeometryTypes::none(dim); // octahedron
        else
            DUNE_THROW(Dune::NotImplemented, "PQ1Bubble scv geometries for dim=" << dim
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
            using Corners = Detail::PQ1Bubble::OverlappingScvfCorners<Dune::GeometryTypes::triangle>;
            return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localOverlappingScvfIdx]);
        }
        else if (type == Dune::GeometryTypes::quadrilateral)
        {
            using Corners = Detail::PQ1Bubble::OverlappingScvfCorners<Dune::GeometryTypes::quadrilateral>;
            return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localOverlappingScvfIdx]);
        }
        else if (type == Dune::GeometryTypes::tetrahedron)
        {
            using Corners = Detail::PQ1Bubble::OverlappingScvfCorners<Dune::GeometryTypes::tetrahedron>;
            return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localOverlappingScvfIdx]);
        }
        else if (type == Dune::GeometryTypes::hexahedron)
        {
            using Corners = Detail::PQ1Bubble::OverlappingScvfCorners<Dune::GeometryTypes::hexahedron>;
            return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(geo_, Corners::keys[localOverlappingScvfIdx]);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "PQ1Bubble scvf geometries for dim=" << dim
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
        return boxHelper_.getBoundaryScvfCorners(localFacetIndex, indexInFacet);
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
        return boxHelper_.numInteriorScvf() + referenceElement(geo_).size(dim);
    }

    //! number of boundary sub control volume faces for face localFacetIndex
    std::size_t numBoundaryScvf(unsigned int localFacetIndex) const
    {
        return referenceElement(geo_).size(localFacetIndex, 1, dim);
    }

    //! number of sub control volumes (number of codim-1 entities)
    std::size_t numScv() const
    {
        return boxHelper_.numScv() + 1;
    }

    //! get scv volume
    Scalar scvVolume(unsigned int localScvIdx, const ScvCornerStorage& p) const
    {
        const auto scvType = getScvGeometryType(localScvIdx);
        if constexpr (dim == 3)
            if (scvType == Dune::GeometryTypes::none(dim))
                return octahedronVolume_(p);

        return Dumux::convexPolytopeVolume<dim>(
            scvType,
            [&](unsigned int i){ return p[i]; }
        );
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
        if (localScvfIndex < numEdges)
            return {
                static_cast<LocalIndexType>(referenceElement(geo_).subEntity(localScvfIndex, dim-1, 0, dim)),
                static_cast<LocalIndexType>(referenceElement(geo_).subEntity(localScvfIndex, dim-1, 1, dim))
            };
        else
            return {
                static_cast<LocalIndexType>(numScv()-1),
                static_cast<LocalIndexType>(localScvfIndex-numEdges)
            };
    }

    std::array<LocalIndexType, 2> getScvPairForBoundaryScvf(unsigned int localFacetIndex, unsigned int localIsScvfIndex) const
    {
        const LocalIndexType insideScvIdx
            = static_cast<LocalIndexType>(referenceElement(geo_).subEntity(localFacetIndex, 1, localIsScvfIndex, dim));
        return { insideScvIdx, insideScvIdx };
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

private:
    Scalar octahedronVolume_(const ScvCornerStorage& p) const
    {
        using std::abs;
        return 1.0/6.0 * (
            abs(Dumux::tripleProduct(p[4]-p[0], p[1]-p[0], p[2]-p[0]))
            + abs(Dumux::tripleProduct(p[5]-p[0], p[1]-p[0], p[2]-p[0]))
        );
    }

    const typename Element::Geometry& geo_; //!< Reference to the element geometry
    Dumux::BoxGeometryHelper<GridView, dim, ScvType, ScvfType> boxHelper_;
};

} // end namespace Dumux

#endif
