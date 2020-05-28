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
 * \ingroup FacetCoupling
 * \brief \copydoc Dumux::BoxFacetCouplingFVElementGeometry
 */
#ifndef DUMUX_FACETCOUPLING_BOX_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_FACETCOUPLING_BOX_FV_ELEMENT_GEOMETRY_HH

#include <algorithm>

#include <dune/geometry/type.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/scvandscvfiterators.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Base class for the element-local finite volume geometry for box models
 *        in the context of models considering coupling of different domains across the
 *        bulk grid facets. This builds up the sub control volumes and sub control volume
 *        faces for an element.
 * \tparam GG the finite volume grid geometry type
 * \tparam enableGridGeometryCache if the grid geometry is cached or not
 */
template<class GG, bool enableGridGeometryCache>
class BoxFacetCouplingFVElementGeometry;

//! specialization in case the FVElementGeometries are stored
template<class GG>
class BoxFacetCouplingFVElementGeometry<GG, true>
{
    using GridView = typename GG::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;
public:
    //! export type of the element
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using GridGeometry = GG;
    //! the maximum number of scvs per element (2^dim for cubes)
    static constexpr std::size_t maxNumElementScvs = (1<<dim);

    //! Constructor
    BoxFacetCouplingFVElementGeometry(const GridGeometry& gridGeometry)
    : gridGeometryPtr_(&gridGeometry)
    {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(LocalIndexType scvIdx) const
    { return gridGeometry().scvs(eIdx_)[scvIdx]; }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(LocalIndexType scvfIdx) const
    { return gridGeometry().scvfs(eIdx_)[scvfIdx]; }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolume>::const_iterator>
    scvs(const BoxFacetCouplingFVElementGeometry& fvGeometry)
    {
        const auto& g = fvGeometry.gridGeometry();
        using Iter = typename std::vector<SubControlVolume>::const_iterator;
        return Dune::IteratorRange<Iter>(g.scvs(fvGeometry.eIdx_).begin(), g.scvs(fvGeometry.eIdx_).end());
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolumeFace>::const_iterator>
    scvfs(const BoxFacetCouplingFVElementGeometry& fvGeometry)
    {
        const auto& g = fvGeometry.gridGeometry();
        using Iter = typename std::vector<SubControlVolumeFace>::const_iterator;
        return Dune::IteratorRange<Iter>(g.scvfs(fvGeometry.eIdx_).begin(), g.scvfs(fvGeometry.eIdx_).end());
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    { return gridGeometry().feCache().get(elemGeometryType_).localBasis(); }

    //! The total number of sub control volumes
    std::size_t numScv() const
    { return gridGeometry().scvs(eIdx_).size(); }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    { return gridGeometry().scvfs(eIdx_).size(); }

    //! this function is for compatibility reasons with cc methods
    //! The box stencil is always element-local so bind and bindElement
    //! are identical.
    void bind(const Element& element)
    {
        this->bindElement(element);
    }

    //! Binding of an element, has to be called before using the fvgeometries
    //! Prepares all the volume variables within the element
    //! For compatibility reasons with the FVGeometry cache being disabled
    void bindElement(const Element& element)
    {
        elemGeometryType_ = element.type();
        eIdx_ = gridGeometry().elementMapper().index(element);
    }

    //! The global finite volume geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return *gridGeometryPtr_; }

private:
    Dune::GeometryType elemGeometryType_;
    const GridGeometry* gridGeometryPtr_;

    GridIndexType eIdx_;
};

//! specialization in case the geometries are not stored grid-wide
template<class GG>
class BoxFacetCouplingFVElementGeometry<GG, false>
{
    using GridView = typename GG::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;

    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;

    using GeometryHelper = BoxGeometryHelper<GridView, dim,
                                             typename GG::SubControlVolume,
                                             typename GG::SubControlVolumeFace>;
public:
    //! export type of the element
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using GridGeometry = GG;
    //! the maximum number of scvs per element (2^dim for cubes)
    static constexpr std::size_t maxNumElementScvs = (1<<dim);

    //! Constructor
    BoxFacetCouplingFVElementGeometry(const GridGeometry& gridGeometry)
    : gridGeometryPtr_(&gridGeometry)
    {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(LocalIndexType scvIdx) const
    { return scvs_[scvIdx]; }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(LocalIndexType scvfIdx) const
    { return scvfs_[scvfIdx]; }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolume>::const_iterator>
    scvs(const BoxFacetCouplingFVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolume>::const_iterator;
        return Dune::IteratorRange<Iter>(fvGeometry.scvs_.begin(), fvGeometry.scvs_.end());
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolumeFace>::const_iterator>
    scvfs(const BoxFacetCouplingFVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolumeFace>::const_iterator;
        return Dune::IteratorRange<Iter>(fvGeometry.scvfs_.begin(), fvGeometry.scvfs_.end());
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    { return gridGeometry().feCache().get(elemGeometryType_).localBasis(); }

    //! The total number of sub control volumes
    std::size_t numScv() const
    { return scvs_.size(); }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    { return scvfs_.size(); }

    //! this function is for compatibility reasons with cc methods
    //! The box stencil is always element-local so bind and bindElement
    //! are identical.
    void bind(const Element& element)
    {
        this->bindElement(element);
    }

    //! Binding of an element, has to be called before using the fvgeometries
    //! Prepares all the volume variables within the element
    //! For compatibility reasons with the FVGeometry cache being disabled
    void bindElement(const Element& element)
    {
        eIdx_ = gridGeometry().elementMapper().index(element);
        makeElementGeometries(element);
    }

    //! The global finite volume geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return *gridGeometryPtr_; }

private:

    void makeElementGeometries(const Element& element)
    {
        auto eIdx = gridGeometry().elementMapper().index(element);

        // get the element geometry
        auto elementGeometry = element.geometry();
        elemGeometryType_ = elementGeometry.type();
        const auto refElement = referenceElement(elementGeometry);

        // get the sub control volume geometries of this element
        GeometryHelper geometryHelper(elementGeometry);

        // construct the sub control volumes
        scvs_.clear();
        scvs_.reserve(elementGeometry.corners());
        using LocalIndexType = typename SubControlVolumeFace::Traits::LocalIndexType;
        for (LocalIndexType scvLocalIdx = 0; scvLocalIdx < elementGeometry.corners(); ++scvLocalIdx)
            scvs_.emplace_back(geometryHelper,
                               scvLocalIdx,
                               eIdx,
                               gridGeometry().vertexMapper().subIndex(element, scvLocalIdx, dim));

        // construct the sub control volume faces
        const auto numInnerScvf = (dim==1) ? 1 : element.subEntities(dim-1);
        scvfs_.clear();
        scvfs_.reserve(numInnerScvf);

        unsigned int scvfLocalIdx = 0;
        for (; scvfLocalIdx < numInnerScvf; ++scvfLocalIdx)
        {
            // find the local scv indices this scvf is connected to
            std::vector<LocalIndexType> localScvIndices({static_cast<LocalIndexType>(refElement.subEntity(scvfLocalIdx, dim-1, 0, dim)),
                                                         static_cast<LocalIndexType>(refElement.subEntity(scvfLocalIdx, dim-1, 1, dim))});

            // create the sub-control volume face
            scvfs_.emplace_back(geometryHelper,
                                element,
                                elementGeometry,
                                scvfLocalIdx,
                                std::move(localScvIndices));
        }

        // construct the sub control volume faces on the domain/interior boundaries
        // skip handled facets (necessary for e.g. Dune::FoamGrid)
        std::vector<unsigned int> handledFacets;
        for (const auto& intersection : intersections(gridGeometry().gridView(), element))
        {
            if (std::count(handledFacets.begin(), handledFacets.end(), intersection.indexInInside()))
                continue;

            handledFacets.push_back(intersection.indexInInside());

            // determine if all corners live on the facet grid
            const auto isGeometry = intersection.geometry();
            const auto numFaceCorners = isGeometry.corners();
            const auto idxInInside = intersection.indexInInside();
            const auto boundary = intersection.boundary();

            std::vector<LocalIndexType> vIndicesLocal(numFaceCorners);
            for (int i = 0; i < numFaceCorners; ++i)
                vIndicesLocal[i] = static_cast<LocalIndexType>(refElement.subEntity(idxInInside, 1, i, dim));

            // if all vertices are living on the facet grid, this is an interiour boundary
            const bool isOnFacet = gridGeometry().isOnInteriorBoundary(element, intersection);

            if (isOnFacet || boundary)
            {
                for (unsigned int isScvfLocalIdx = 0; isScvfLocalIdx < isGeometry.corners(); ++isScvfLocalIdx)
                {
                    // find the inside scv this scvf is belonging to (localIdx = element local vertex index)
                    std::vector<LocalIndexType> localScvIndices = {vIndicesLocal[isScvfLocalIdx], vIndicesLocal[isScvfLocalIdx]};

                    // create the sub-control volume face
                    scvfs_.emplace_back(geometryHelper,
                                        intersection,
                                        isGeometry,
                                        isScvfLocalIdx,
                                        scvfLocalIdx,
                                        std::move(localScvIndices),
                                        boundary,
                                        isOnFacet);

                    // increment local counter
                    scvfLocalIdx++;
                }
            }
        }
    }

    //! The bound element
    Dune::GeometryType elemGeometryType_;
    GridIndexType eIdx_;

    //! The global geometry this is a restriction of
    const GridGeometry* gridGeometryPtr_;

    //! vectors to store the geometries locally after binding an element
    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;
};

} // end namespace Dumux

#endif
