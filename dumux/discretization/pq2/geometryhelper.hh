// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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

    using BoxHelper = Dumux::BoxGeometryHelper<GridView, dim, ScvType, ScvfType>;
public:

    HybridPQ2GeometryHelper(const typename Element::Geometry& geometry)
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

    Dune::GeometryType getBoundaryScvfGeometryType(unsigned int localScvfIdx) const
    {
        return Dune::GeometryTypes::cube(dim-1);
    }

    template<int d = dimWorld, std::enable_if_t<(d==3), int> = 0>
    GlobalPosition normal(const ScvfCornerStorage& p, const std::array<LocalIndexType, 2>& scvPair)
    {
        auto normal = Dumux::crossProduct(p[1]-p[0], p[2]-p[0]);
        normal /= normal.two_norm();

        GlobalPosition v = geo_.corner(scvPair[1]) - geo_.corner(scvPair[0]);

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

        GlobalPosition v = geo_.corner(scvPair[1]) - geo_.corner(scvPair[0]);

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

    //! Local dof index related to a localDof, with index ilocalDofIdx, on an intersection with index iIdx
    template<class LocalKey>
    static auto localDofOnIntersection(Dune::GeometryType type, unsigned int iIdx, const LocalKey& localKey)
    {
        const auto& refElement = Dune::referenceElement<Scalar, dim>(type);

        const auto numEntitiesIntersection = refElement.size(iIdx, 1, localKey.codim());
        for(std::size_t idx=0; idx < numEntitiesIntersection; idx++)
            if(localKey.subEntity() == refElement.subEntity(iIdx, 1, idx, localKey.codim()))
                return true;

        return false;
    }

    template<class DofMapper, class LocalKey>
    static auto dofIndex(const DofMapper& dofMapper, const Element& element, const LocalKey& localKey)
    {
        // All dofs are directly related to grid entities, i.e localKey.index() is always zero
        return dofMapper.subIndex(element, localKey.subEntity(), localKey.codim());
    }

    template<class Geometry, class LocalKey>
    static GlobalPosition dofPosition(const Geometry& geo, const LocalKey& localKey)
    {
        if(localKey.codim() == dim)
            return geo.corner(localKey.subEntity());
        else if(localKey.codim() == 0) // should only be called for cubes
            return geo.center();
        else
            return geo.global(localDofPosition(geo.type(), localKey));
    }

    template<class LocalKey>
    GlobalPosition dofPosition(const LocalKey& localKey) const
    {
        return dofPosition(geo_, localKey);
    }

    //! local dof position
    template<class LocalKey>
    static Element::Geometry::LocalCoordinate localDofPosition(Dune::GeometryType type, const LocalKey& localKey)
    {
        return Dune::referenceElement<Scalar, dim>(type).position(localKey.subEntity(), localKey.codim());
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

    // For hybrid dofs we don't construct scvfs
    bool isOverlappingScvf(unsigned int localScvfIndex) const
    { return false; }

    bool isOverlappingBoundaryScvf(unsigned int localFacetIndex) const
    { return false; }

    // For hybrid dofs we don't construct scvs
    bool isOverlappingScv(unsigned int localScvIndex) const
    { return false; }

private:
    const typename Element::Geometry& geo_; //!< Reference to the element geometry
    BoxHelper boxHelper_;
};

} // end namespace Dumux

#endif
