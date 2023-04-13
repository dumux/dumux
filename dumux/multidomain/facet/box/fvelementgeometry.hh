// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetCoupling
 * \brief \copydoc Dumux::BoxFacetCouplingFVElementGeometry
 */
#ifndef DUMUX_FACETCOUPLING_BOX_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_FACETCOUPLING_BOX_FV_ELEMENT_GEOMETRY_HH

#include <algorithm>
#include <optional>

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
    { return gridGeometry().feCache().get(element_->type()).localBasis(); }

    //! The total number of sub control volumes
    std::size_t numScv() const
    { return gridGeometry().scvs(eIdx_).size(); }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    { return gridGeometry().scvfs(eIdx_).size(); }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    BoxFacetCouplingFVElementGeometry bind(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! this function is for compatibility reasons with cc methods
    //! The box stencil is always element-local so bind and bindElement
    //! are identical.
    void bind(const Element& element) &
    { this->bindElement(element); }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bindElement(element);`
     */
    BoxFacetCouplingFVElementGeometry bindElement(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! Binding of an element, has to be called before using the fvgeometries
    //! Prepares all the volume variables within the element
    //! For compatibility reasons with the FVGeometry cache being disabled
    void bindElement(const Element& element) &
    {
        element_ = element;
        eIdx_ = gridGeometry().elementMapper().index(element);
    }

    //! The global finite volume geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return *gridGeometryPtr_; }

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return static_cast<bool>(element_); }

    //! The bound element
    const Element& element() const
    { return *element_; }

    // suppress warnings due to current implementation
    // these interfaces should be used!
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"

    //! Create the geometry of a given sub control volume
    typename SubControlVolume::Traits::Geometry geometry(const SubControlVolume& scv) const
    { return scv.geometry(); }

    //! Create the geometry of a given sub control volume face
    typename SubControlVolumeFace::Traits::Geometry geometry(const SubControlVolumeFace& scvf) const
    { return scvf.geometry(); }

    #pragma GCC diagnostic pop

private:
    const GridGeometry* gridGeometryPtr_;

    GridIndexType eIdx_;
    std::optional<Element> element_;
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

    using GeometryHelper = typename GG::GeometryHelper;
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
    { return gridGeometry().feCache().get(element_->type()).localBasis(); }

    //! The total number of sub control volumes
    std::size_t numScv() const
    { return scvs_.size(); }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    { return scvfs_.size(); }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    BoxFacetCouplingFVElementGeometry bind(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! this function is for compatibility reasons with cc methods
    //! The box stencil is always element-local so bind and bindElement
    //! are identical.
    void bind(const Element& element) &
    { this->bindElement(element); }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bindElement(element);`
     */
    BoxFacetCouplingFVElementGeometry bindElement(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! Binding of an element, has to be called before using the fvgeometries
    //! Prepares all the volume variables within the element
    //! For compatibility reasons with the FVGeometry cache being disabled
    void bindElement(const Element& element) &
    {
        element_ = element;
        eIdx_ = gridGeometry().elementMapper().index(element);
        makeElementGeometries_();
    }

    //! The global finite volume geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return *gridGeometryPtr_; }

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return static_cast<bool>(element_); }

    //! The bound element
    const Element& element() const
    { return *element_; }

    // suppress warnings due to current implementation
    // these interfaces should be used!
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"

    //! Create the geometry of a given sub control volume
    typename SubControlVolume::Traits::Geometry geometry(const SubControlVolume& scv) const
    { return scv.geometry(); }

    //! Create the geometry of a given sub control volume face
    typename SubControlVolumeFace::Traits::Geometry geometry(const SubControlVolumeFace& scvf) const
    { return scvf.geometry(); }

    #pragma GCC diagnostic pop

private:

    void makeElementGeometries_()
    {
        // get the element geometry
        const auto& element = *element_;
        const auto elementGeometry = element.geometry();
        const auto refElement = referenceElement(elementGeometry);

        // get the sub control volume geometries of this element
        GeometryHelper geometryHelper(elementGeometry);

        // construct the sub control volumes
        scvs_.clear();
        scvs_.reserve(elementGeometry.corners());
        using LocalIndexType = typename SubControlVolumeFace::Traits::LocalIndexType;
        for (LocalIndexType scvLocalIdx = 0; scvLocalIdx < elementGeometry.corners(); ++scvLocalIdx)
            scvs_.emplace_back(geometryHelper.getScvCorners(scvLocalIdx),
                               scvLocalIdx,
                               eIdx_,
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

    //! The global geometry this is a restriction of
    const GridGeometry* gridGeometryPtr_;

    //! The bound element
    GridIndexType eIdx_;
    std::optional<Element> element_;

    //! vectors to store the geometries locally after binding an element
    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;
};

} // end namespace Dumux

#endif
