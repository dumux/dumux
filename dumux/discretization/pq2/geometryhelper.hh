// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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

/*!
 * \ingroup PQ2Discretization
 * \brief A class to create sub control volume and sub control volume face geometries per element
 */
template <class GridView, class ScvType, class ScvfType>
class HybridPQ2GeometryHelper
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

    HybridPQ2GeometryHelper(const typename Element::Geometry& geometry)
    : geo_(geometry)
    , boxHelper_(geometry)
    {}

    //! Create a vector with the scv corners
    ScvCornerStorage getScvCorners(unsigned int localScvIdx) const
    {
        // proceed according to number of corners of the element
        const auto numBoxScv = boxHelper_.numScv();
        // reuse box geometry helper for the corner scvs
        if (localScvIdx < numBoxScv)
            return boxHelper_.getScvCorners(localScvIdx);

        DUNE_THROW(Dune::NotImplemented, "PQ2 scv corners call for hybrid dofs");
    }

    Dune::GeometryType getScvGeometryType(unsigned int localScvIdx) const
    {
        // proceed according to number of corners of the element
        const auto numBoxScv = boxHelper_.numScv();

        if (localScvIdx < numBoxScv)
            return Dune::GeometryTypes::cube(dim);

        DUNE_THROW(Dune::NotImplemented, "PQ2 scv geometry call for hybrid dofs");
    }

    //! Create a vector with the corners of sub control volume faces
    ScvfCornerStorage getScvfCorners(unsigned int localScvfIdx) const
    {
        // proceed according to number of corners
        const auto numBoxScvf = boxHelper_.numInteriorScvf();
        // reuse box geometry helper for the corner scvs
        if (localScvfIdx < numBoxScvf)
            return boxHelper_.getScvfCorners(localScvfIdx);

        DUNE_THROW(Dune::NotImplemented, "PQ2 scvf corners call for hybrid dofs");
    }

    Dune::GeometryType getInteriorScvfGeometryType(unsigned int localScvfIdx) const
    {
        const auto numBoxScvf = boxHelper_.numInteriorScvf();
        if (localScvfIdx < numBoxScvf)
            return Dune::GeometryTypes::cube(dim-1);

        DUNE_THROW(Dune::NotImplemented, "PQ2 interior scvf geometry type call for hybrid dofs");
    }

    //! Create the sub control volume face geometries on the boundary
    ScvfCornerStorage getBoundaryScvfCorners(unsigned int localFacetIndex,
                                             unsigned int indexInFacet) const
    {
        return boxHelper_.getBoundaryScvfCorners(localFacetIndex, indexInFacet);
    }

    Dune::GeometryType getBoundaryScvfGeometryType(unsigned int localFacetIndex,
                                                   unsigned int indexInFacet) const
    {
        const auto numBoxScvf = Dune::referenceElement(geo_).size(localFacetIndex, 1, dim);
        if (indexInFacet < numBoxScvf)
            return Dune::GeometryTypes::cube(dim-1);

        DUNE_THROW(Dune::NotImplemented, "PQ2 boundary scvf type for hybrid faces: " << "indexInFacet: " << indexInFacet);
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
        return boxHelper_.numInteriorScvf();
    }

    //! number of boundary sub control volume faces for face localFacetIndex
    std::size_t numBoundaryScvf(unsigned int localFacetIndex) const
    {
        return referenceElement(geo_).size(localFacetIndex, 1, dim);
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
    std::size_t numElementDofs() const
    {
        return boxHelper_.numScv() + referenceElement(geo_).size(dim-1);
    }

    //! number of hybrid dofs
    static std::size_t numHybridDofs(Dune::GeometryType type)
    {
        const auto ref = referenceElement<double,dim>(type);
        return ref.size(dim-1);
    }

    template<class DofMapper>
    static auto dofIndex(const DofMapper& dofMapper, const Element& element, unsigned int localScvIdx)
    {
        const auto numBoxScv = Dune::referenceElement<Scalar, dim>(element.type()).size(dim);
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

        DUNE_THROW(Dune::NotImplemented, "PQ2 scv pair call for hybrid dofs");
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

        DUNE_THROW(Dune::NotImplemented, "PQ2 scv boundary pair call for hybrid dofs");
    }

    bool isOverlappingBoundaryScvf(unsigned int localFacetIndex, unsigned int localIsScvfIndex) const
    {
        const auto numBoxScvf = referenceElement(geo_).size(localFacetIndex, 1, dim);
        if (localIsScvfIndex < numBoxScvf)
            return false;
        else
            return true;
    }

    // For hybrid scheme we only have non-overlapping entities
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
    LocalIndexType localKeyToLocalDofIndex(const LocalKey& localKey) const
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
