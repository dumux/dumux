// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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
#include <dumux/geometry/center.hh>

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

    using BoxHelper = Dumux::BoxGeometryHelper<GridView, dim, ScvType, ScvfType>;
public:

    PQ1BubbleGeometryHelper(const typename Element::Geometry& geometry)
    : geo_(geometry)
    , boxHelper_(geometry)
    {}

    //! Create a vector with the scv corners
    ScvCornerStorage getScvCorners(unsigned int localScvIdx) const
    {
        return getScvCorners(geo_.type(), [&](const auto& local){ return geo_.global(local); },  localScvIdx);
    }

    //! Create a vector with the scv corners
    template<class Transformation>
    static ScvCornerStorage getScvCorners(Dune::GeometryType type, Transformation&& trans, unsigned int localScvIdx)
    {
        // proceed according to number of corners of the element
        const auto& ref = Dune::referenceElement<Scalar, dim>(type);
        const auto numBoxScv = ref.size(dim);
        // reuse box geometry helper for the corner scvs
        if (localScvIdx < numBoxScv)
            return BoxHelper::getScvCorners(type, trans, localScvIdx);

        const auto localOverlappingScvIdx = localScvIdx-numBoxScv;
        if (type == Dune::GeometryTypes::triangle)
        {
            using Corners = Detail::PQ1Bubble::OverlappingScvCorners<Dune::GeometryTypes::triangle>;
            return Detail::Box::keyToCornerStorage<ScvCornerStorage>(ref, trans, Corners::keys[localOverlappingScvIdx]);
        }
        else if (type == Dune::GeometryTypes::quadrilateral)
        {
            using Corners = Detail::PQ1Bubble::OverlappingScvCorners<Dune::GeometryTypes::quadrilateral>;
            return Detail::Box::keyToCornerStorage<ScvCornerStorage>(ref, trans, Corners::keys[localOverlappingScvIdx]);
        }
        else if (type == Dune::GeometryTypes::tetrahedron)
        {
            using Corners = Detail::PQ1Bubble::OverlappingScvCorners<Dune::GeometryTypes::tetrahedron>;
            return Detail::Box::keyToCornerStorage<ScvCornerStorage>(ref, trans, Corners::keys[localOverlappingScvIdx]);
        }
        else if (type == Dune::GeometryTypes::hexahedron)
        {
            using Corners = Detail::PQ1Bubble::OverlappingScvCorners<Dune::GeometryTypes::hexahedron>;
            return Detail::Box::keyToCornerStorage<ScvCornerStorage>(ref, trans, Corners::keys[localOverlappingScvIdx]);
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
        return getScvfCorners(geo_.type(), [&](const auto& local){ return geo_.global(local); },  localScvfIdx);
    }

    //! Create a vector with the corners of sub control volume faces
    template<class Transformation>
    static ScvfCornerStorage getScvfCorners(Dune::GeometryType type, Transformation&& trans, unsigned int localScvfIdx)
    {
        // proceed according to number of corners
        const auto& ref = Dune::referenceElement<Scalar, dim>(type);
        const auto numBoxScvf = ref.size(dim-1);
        // reuse box geometry helper for the corner scvs
        if (localScvfIdx < numBoxScvf)
            return BoxHelper::getScvfCorners(type, trans, localScvfIdx);

        const auto localOverlappingScvfIdx = localScvfIdx-numBoxScvf;
        if (type == Dune::GeometryTypes::triangle)
        {
            using Corners = Detail::PQ1Bubble::OverlappingScvfCorners<Dune::GeometryTypes::triangle>;
            return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(ref, trans, Corners::keys[localOverlappingScvfIdx]);
        }
        else if (type == Dune::GeometryTypes::quadrilateral)
        {
            using Corners = Detail::PQ1Bubble::OverlappingScvfCorners<Dune::GeometryTypes::quadrilateral>;
            return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(ref, trans, Corners::keys[localOverlappingScvfIdx]);
        }
        else if (type == Dune::GeometryTypes::tetrahedron)
        {
            using Corners = Detail::PQ1Bubble::OverlappingScvfCorners<Dune::GeometryTypes::tetrahedron>;
            return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(ref, trans, Corners::keys[localOverlappingScvfIdx]);
        }
        else if (type == Dune::GeometryTypes::hexahedron)
        {
            using Corners = Detail::PQ1Bubble::OverlappingScvfCorners<Dune::GeometryTypes::hexahedron>;
            return Detail::Box::keyToCornerStorage<ScvfCornerStorage>(ref, trans, Corners::keys[localOverlappingScvfIdx]);
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

    Dune::GeometryType getBoundaryScvfGeometryType(unsigned int localScvfIdx) const
    {
        return Dune::GeometryTypes::cube(dim-1);
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
    static auto numInteriorScvf(Dune::GeometryType type)
    {
        return BoxHelper::numInteriorScvf(type) + Dune::referenceElement<Scalar, dim>(type).size(dim);
    }

    //! number of boundary sub control volume faces for face localFacetIndex
    static auto numBoundaryScvf(Dune::GeometryType type, unsigned int localFacetIndex)
    {
        return Dune::referenceElement<Scalar, dim>(type).size(localFacetIndex, 1, dim);
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

    //! number of element dofs
    static std::size_t numElementDofs(Dune::GeometryType type)
    {
        return Dune::referenceElement<Scalar, dim>(type).size(dim) + 1;
    }

    //! number of hybrid dofs
    static std::size_t numNonCVLocalDofs(Dune::GeometryType type)
    {
        return 0;
    }

    //! Number of local dofs related to an intersection with index iIdx
    static auto numLocalDofsIntersection(Dune::GeometryType type, unsigned int iIdx)
    {
        return Dune::referenceElement<Scalar, dim>(type).size(iIdx, 1, dim);
    }

    //! Local dof index related to a localDof, with index ilocalDofIdx, on an intersection with index iIdx
    static auto localDofIndexIntersection(Dune::GeometryType type, unsigned int iIdx, unsigned int ilocalDofIdx)
    {
        return Dune::referenceElement<Scalar, dim>(type).subEntity(iIdx, 1, ilocalDofIdx, dim);
    }

    template<class DofMapper>
    static auto dofIndex(const DofMapper& dofMapper, const Element& element, unsigned int localDofIdx)
    {
        if (localDofIdx < numElementDofs(element.type())-1)
            return dofMapper.subIndex(element, localDofIdx, dim);
        else
            return dofMapper.index(element);
    }

    static GlobalPosition dofPosition(const Element& element, unsigned int localDofIdx)
    {
        if (localDofIdx < numElementDofs(element.geometry().type())-1)
            return element.geometry().corner(localDofIdx);
        else
            return element.geometry().center();
    }

    GlobalPosition dofPosition(unsigned int localDofIdx) const
    {
        if (localDofIdx < numElementDofs(geo_.type())-1)
            return geo_.corner(localDofIdx);
        else
            return geo_.center();
    }

    //! local dof position
    template<class LocalKey>
    static Element::Geometry::LocalCoordinate localDofPosition(Dune::GeometryType type, const LocalKey& localKey)
    {
        return Dune::referenceElement<Scalar, dim>(type).position(localKey.subEntity(), localKey.codim());
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
                static_cast<LocalIndexType>(numElementDofs(geo_.type())-1),
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

    bool isOverlappingBoundaryScvf(unsigned int localFacetIndex) const
    {
        return false;
    }

    bool isOverlappingScv(unsigned int localScvIndex) const
    {
        if (localScvIndex < boxHelper_.numScv())
            return false;
        else
            return true;
    }

    //! local scvf center
    static Element::Geometry::LocalCoordinate localScvfCenter(Dune::GeometryType type, unsigned int localScvfIdx)
    {
        return Dumux::center(getScvfCorners(type, [&](const auto& local){ return local; }, localScvfIdx));
    }

    //! local boundary scvf center
    static Element::Geometry::LocalCoordinate localBoundaryScvfCenter(Dune::GeometryType type, unsigned int localFacetIndex, unsigned int indexInFace)
    {
        return Dumux::center(BoxHelper::getBoundaryScvfCorners(type, [&](const auto& local){ return local; }, localFacetIndex, indexInFace));
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
    BoxHelper boxHelper_;
};

template <class GridView, class ScvType, class ScvfType>
class HybridPQ1BubbleGeometryHelper
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

    using BoxHelper = Dumux::BoxGeometryHelper<GridView, dim, ScvType, ScvfType>;
public:

    HybridPQ1BubbleGeometryHelper(const typename Element::Geometry& geometry)
    : geo_(geometry)
    , boxHelper_(geometry)
    {}

    //! Create a vector with the scv corners
    ScvCornerStorage getScvCorners(unsigned int localScvIdx) const
    {
        return getScvCorners(geo_.type(), [&](const auto& local){ return geo_.global(local); },  localScvIdx);
    }

    //! Create a vector with the scv corners
    template<class Transformation>
    static ScvCornerStorage getScvCorners(Dune::GeometryType type, Transformation&& trans, unsigned int localScvIdx)
    {
        // proceed according to number of corners of the element
        const auto& ref = Dune::referenceElement<Scalar, dim>(type);
        const auto numBoxScv = ref.size(dim);
        // reuse box geometry helper for the corner scvs
        if (localScvIdx < numBoxScv)
            return BoxHelper::getScvCorners(type, trans, localScvIdx);

        DUNE_THROW(Dune::NotImplemented, "PQ1Bubble scv corners call for hybrid dofs");
    }

    Dune::GeometryType getScvGeometryType(unsigned int localScvIdx) const
    {
        // proceed according to number of corners of the element
        const auto numBoxScv = boxHelper_.numScv();

        if (localScvIdx < numBoxScv)
            return Dune::GeometryTypes::cube(dim);

        DUNE_THROW(Dune::NotImplemented, "PQ1Bubble scv geometry call for hybrid dofs");
    }

    //! Create a vector with the corners of sub control volume faces
    ScvfCornerStorage getScvfCorners(unsigned int localScvfIdx) const
    {
        return getScvfCorners(geo_.type(), [&](const auto& local){ return geo_.global(local); },  localScvfIdx);
    }

    //! Create a vector with the corners of sub control volume faces
    template<class Transformation>
    static ScvfCornerStorage getScvfCorners(Dune::GeometryType type, Transformation&& trans, unsigned int localScvfIdx)
    {
        // proceed according to number of corners
        const auto& ref = Dune::referenceElement<Scalar, dim>(type);
        const auto numBoxScvf = ref.size(dim-1);
        // reuse box geometry helper for scvfs
        if (localScvfIdx < numBoxScvf)
            return BoxHelper::getScvfCorners(type, trans, localScvfIdx);

        DUNE_THROW(Dune::NotImplemented, "PQ1Bubble scvf corners call for hybrid dofs");
    }

    Dune::GeometryType getInteriorScvfGeometryType(unsigned int localScvfIdx) const
    {
        const auto numBoxScvf = boxHelper_.numInteriorScvf();
        if (localScvfIdx < numBoxScvf)
            return Dune::GeometryTypes::cube(dim-1);

        DUNE_THROW(Dune::NotImplemented, "PQ1Bubble interior scvf geometry type call for hybrid dofs");
    }

    //! Create the sub control volume face geometries on the boundary
    ScvfCornerStorage getBoundaryScvfCorners(unsigned int localFacetIndex,
                                             unsigned int indexInFacet) const
    {
        return boxHelper_.getBoundaryScvfCorners(localFacetIndex, indexInFacet);
    }

    Dune::GeometryType getBoundaryScvfGeometryType(unsigned int localScvfIdx) const
    {
        return Dune::GeometryTypes::cube(dim-1);
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
    static auto numInteriorScvf(Dune::GeometryType type)
    {
        return BoxHelper::numInteriorScvf(type);
    }

    //! number of boundary sub control volume faces for face localFacetIndex
    static auto numBoundaryScvf(Dune::GeometryType type, unsigned int localFacetIndex)
    {
        return Dune::referenceElement<Scalar, dim>(type).size(localFacetIndex, 1, dim);
    }

    //! number of sub control volumes (number of codim-1 entities)
    std::size_t numScv() const
    {
        return boxHelper_.numScv();
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

    //! number of element dofs
    static std::size_t numElementDofs(Dune::GeometryType type)
    {
        return Dune::referenceElement<Scalar, dim>(type).size(dim) + 1;
    }

    //! number of hybrid dofs
    static std::size_t numNonCVLocalDofs(Dune::GeometryType type)
    {
        return 1;
    }

    //! Number of local dofs related to an intersection with index iIdx
    static auto numLocalDofsIntersection(Dune::GeometryType type, unsigned int iIdx)
    {
        return Dune::referenceElement<Scalar, dim>(type).size(iIdx, 1, dim);
    }

    //! Local dof index related to a localDof, with index ilocalDofIdx, on an intersection with index iIdx
    static auto localDofIndexIntersection(Dune::GeometryType type, unsigned int iIdx, unsigned int ilocalDofIdx)
    {
        return Dune::referenceElement<Scalar, dim>(type).subEntity(iIdx, 1, ilocalDofIdx, dim);
    }

    template<class DofMapper>
    static auto dofIndex(const DofMapper& dofMapper, const Element& element, unsigned int localDofIdx)
    {
        if (localDofIdx < numElementDofs(element.type())-1)
            return dofMapper.subIndex(element, localDofIdx, dim);
        else
            return dofMapper.index(element);
    }

    static GlobalPosition dofPosition(const Element& element, unsigned int localDofIdx)
    {
        if (localDofIdx < numElementDofs(element.geometry().type())-1)
            return element.geometry().corner(localDofIdx);
        else
            return element.geometry().center();
    }

    GlobalPosition dofPosition(unsigned int localDofIdx) const
    {
        if (localDofIdx < numElementDofs(geo_.type())-1)
            return geo_.corner(localDofIdx);
        else
            return geo_.center();
    }

    //! local dof position
    template<class LocalKey>
    static Element::Geometry::LocalCoordinate localDofPosition(Dune::GeometryType type, const LocalKey& localKey)
    {
        return Dune::referenceElement<Scalar, dim>(type).position(localKey.subEntity(), localKey.codim());
    }

    std::array<LocalIndexType, 2> getScvPairForScvf(unsigned int localScvfIndex) const
    {
        const auto numEdges = referenceElement(geo_).size(dim-1);
        if (localScvfIndex < numEdges)
            return {
                static_cast<LocalIndexType>(referenceElement(geo_).subEntity(localScvfIndex, dim-1, 0, dim)),
                static_cast<LocalIndexType>(referenceElement(geo_).subEntity(localScvfIndex, dim-1, 1, dim))
            };

        DUNE_THROW(Dune::NotImplemented, "PQ1Bubble scv pair call for hybrid dofs");
    }

    std::array<LocalIndexType, 2> getScvPairForBoundaryScvf(unsigned int localFacetIndex, unsigned int localIsScvfIndex) const
    {
        const LocalIndexType insideScvIdx
            = static_cast<LocalIndexType>(referenceElement(geo_).subEntity(localFacetIndex, 1, localIsScvfIndex, dim));
        return { insideScvIdx, insideScvIdx };
    }

    bool isOverlappingScvf(unsigned int localScvfIndex) const
    { return false; }

    bool isOverlappingBoundaryScvf(unsigned int localFacetIndex) const
    { return false; }

    bool isOverlappingScv(unsigned int localScvIndex) const
    { return false; }

    //! local scvf center
    static Element::Geometry::LocalCoordinate localScvfCenter(Dune::GeometryType type, unsigned int localScvfIdx)
    {
        return Dumux::center(getScvfCorners(type, [&](const auto& local){ return local; }, localScvfIdx));
    }

    //! local boundary scvf center
    static Element::Geometry::LocalCoordinate localBoundaryScvfCenter(Dune::GeometryType type, unsigned int localFacetIndex, unsigned int indexInFace)
    {
        return Dumux::center(BoxHelper::getBoundaryScvfCorners(type, [&](const auto& local){ return local; }, localFacetIndex, indexInFace));
    }

private:
    const typename Element::Geometry& geo_; //!< Reference to the element geometry
    Dumux::BoxGeometryHelper<GridView, dim, ScvType, ScvfType> boxHelper_;
};

} // end namespace Dumux

#endif
