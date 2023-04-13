// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxDiscretization
 * \brief Base class for the local finite volume geometry for box models
 *        This builds up the sub control volumes and sub control volume faces
 *        for an element.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_BOX_FV_ELEMENT_GEOMETRY_HH

#include <optional>
#include <utility>
#include <unordered_map>
#include <array>
#include <vector>

#include <dune/geometry/type.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/scvandscvfiterators.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup BoxDiscretization
 * \brief Base class for the finite volume geometry vector for box models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 * \tparam GG the finite volume grid geometry type
 * \tparam enableGridGeometryCache if the grid geometry is cached or not
 */
template<class GG, bool enableGridGeometryCache>
class BoxFVElementGeometry;

//! specialization in case the FVElementGeometries are stored
template<class GG>
class BoxFVElementGeometry<GG, true>
{
    using GridView = typename GG::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;
    using GGCache = typename GG::Cache;
    using GeometryHelper = typename GGCache::GeometryHelper;
public:
    //! export the element type
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using GridGeometry = GG;
    //! the maximum number of scvs per element (2^dim for cubes)
    static constexpr std::size_t maxNumElementScvs = (1<<dim);

    /*!
     * \brief Constructor
     * \note Never use this directly and always construct this class via `localView(gridGeometry)`
     */
    BoxFVElementGeometry(const GGCache& ggCache)
    : ggCache_(&ggCache) {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(LocalIndexType scvIdx) const
    {
        return ggCache_->scvs(eIdx_)[scvIdx];
    }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(LocalIndexType scvfIdx) const
    {
        return ggCache_->scvfs(eIdx_)[scvfIdx];
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolume>::const_iterator>
    scvs(const BoxFVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolume>::const_iterator;
        const auto& s = fvGeometry.ggCache_->scvs(fvGeometry.eIdx_);
        return Dune::IteratorRange<Iter>(s.begin(), s.end());
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolumeFace>::const_iterator>
    scvfs(const BoxFVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolumeFace>::const_iterator;
        const auto& s = fvGeometry.ggCache_->scvfs(fvGeometry.eIdx_);
        return Dune::IteratorRange<Iter>(s.begin(), s.end());
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    {
        return gridGeometry().feCache().get(element_->type()).localBasis();
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return ggCache_->scvs(eIdx_).size();
    }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    {
        return ggCache_->scvfs(eIdx_).size();
    }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    BoxFVElementGeometry bind(const Element& element) &&
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
    BoxFVElementGeometry bindElement(const Element& element) &&
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
        // cache element index
        eIdx_ = gridGeometry().elementMapper().index(element);
    }

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return static_cast<bool>(element_); }

    //! The bound element
    const Element& element() const
    { return *element_; }

    //! The bound element's index in the grid view
    GridIndexType elementIndex() const
    { return eIdx_; }

    //! The grid geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return ggCache_->gridGeometry(); }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return ggCache_->hasBoundaryScvf(eIdx_); }

    //! Geometry of a sub control volume
    typename SubControlVolume::Traits::Geometry geometry(const SubControlVolume& scv) const
    {
        assert(isBound());
        const auto geo = element().geometry();
        return { Dune::GeometryTypes::cube(dim), GeometryHelper(geo).getScvCorners(scv.indexInElement()) };
    }

    //! Geometry of a sub control volume face
    typename SubControlVolumeFace::Traits::Geometry geometry(const SubControlVolumeFace& scvf) const
    {
        assert(isBound());
        const auto geo = element().geometry();
        const GeometryHelper geometryHelper(geo);
        if (scvf.boundary())
        {
            const auto localBoundaryIndex = scvf.index() - geometryHelper.numInteriorScvf();
            const auto& key = ggCache_->scvfBoundaryGeometryKeys(eIdx_)[localBoundaryIndex];
            return { Dune::GeometryTypes::cube(dim-1), geometryHelper.getBoundaryScvfCorners(key[0], key[1]) };
        }
        else
            return { Dune::GeometryTypes::cube(dim-1), geometryHelper.getScvfCorners(scvf.index()) };
    }

private:
    const GGCache* ggCache_;
    GridIndexType eIdx_;

    std::optional<Element> element_;
};

//! specialization in case the FVElementGeometries are not stored
template<class GG>
class BoxFVElementGeometry<GG, false>
{
    using GridView = typename GG::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;
    using GGCache = typename GG::Cache;
    using GeometryHelper = typename GGCache::GeometryHelper;
public:
    //! export the element type
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using GridGeometry = GG;
    //! the maximum number of scvs per element (2^dim for cubes)
    static constexpr std::size_t maxNumElementScvs = (1<<dim);

    /*!
     * \brief Constructor
     * \note Never use this directly and always construct this class via `localView(gridGeometry)`
     */
    BoxFVElementGeometry(const GGCache& ggCache)
    : ggCache_(&ggCache) {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(LocalIndexType scvIdx) const
    {
        return scvs_[scvIdx];
    }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(LocalIndexType scvfIdx) const
    {
        return scvfs_[scvfIdx];
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolume>::const_iterator>
    scvs(const BoxFVElementGeometry& fvGeometry)
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
    scvfs(const BoxFVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolumeFace>::const_iterator;
        return Dune::IteratorRange<Iter>(fvGeometry.scvfs_.begin(), fvGeometry.scvfs_.end());
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    {
        return gridGeometry().feCache().get(element_->type()).localBasis();
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return scvs_.size();
    }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    {
        return scvfs_.size();
    }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    BoxFVElementGeometry bind(const Element& element) &&
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
    BoxFVElementGeometry bindElement(const Element& element) &&
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

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return static_cast<bool>(element_); }

    //! The bound element
    const Element& element() const
    { return *element_; }

    //! The bound element's index in the grid view
    GridIndexType elementIndex() const
    { return eIdx_; }

    //! The grid geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return ggCache_->gridGeometry(); }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return hasBoundaryScvf_; }

    //! Geometry of a sub control volume
    typename SubControlVolume::Traits::Geometry geometry(const SubControlVolume& scv) const
    {
        assert(isBound());
        const auto geo = element().geometry();
        return { Dune::GeometryTypes::cube(dim), GeometryHelper(geo).getScvCorners(scv.indexInElement()) };
    }

    //! Geometry of a sub control volume face
    typename SubControlVolumeFace::Traits::Geometry geometry(const SubControlVolumeFace& scvf) const
    {
        assert(isBound());
        const auto geo = element().geometry();
        const GeometryHelper geometryHelper(geo);
        if (scvf.boundary())
        {
            const auto localBoundaryIndex = scvf.index() - geometryHelper.numInteriorScvf();
            const auto& key = scvfBoundaryGeometryKeys_[localBoundaryIndex];
            return { Dune::GeometryTypes::cube(dim-1), geometryHelper.getBoundaryScvfCorners(key[0], key[1]) };
        }
        else
            return { Dune::GeometryTypes::cube(dim-1), geometryHelper.getScvfCorners(scvf.index()) };
    }

private:
    void makeElementGeometries_()
    {
        hasBoundaryScvf_ = false;

        // get the element geometry
        const auto& element = *element_;
        const auto elementGeometry = element.geometry();
        const auto refElement = referenceElement(elementGeometry);

        // get the sub control volume geometries of this element
        GeometryHelper geometryHelper(elementGeometry);

        // construct the sub control volumes
        scvs_.resize(elementGeometry.corners());
        for (LocalIndexType scvLocalIdx = 0; scvLocalIdx < elementGeometry.corners(); ++scvLocalIdx)
        {
            // get associated dof index
            const auto dofIdxGlobal = gridGeometry().vertexMapper().subIndex(element, scvLocalIdx, dim);

            // add scv to the local container
            scvs_[scvLocalIdx] = SubControlVolume(
                geometryHelper.getScvCorners(scvLocalIdx),
                scvLocalIdx,
                eIdx_,
                dofIdxGlobal
            );
        }

        // construct the sub control volume faces
        const auto numInnerScvf = geometryHelper.numInteriorScvf();
        scvfs_.resize(numInnerScvf);
        scvfBoundaryGeometryKeys_.clear();

        LocalIndexType scvfLocalIdx = 0;
        for (; scvfLocalIdx < numInnerScvf; ++scvfLocalIdx)
        {
            // find the local scv indices this scvf is connected to
            std::vector<LocalIndexType> localScvIndices({static_cast<LocalIndexType>(refElement.subEntity(scvfLocalIdx, dim-1, 0, dim)),
                                                         static_cast<LocalIndexType>(refElement.subEntity(scvfLocalIdx, dim-1, 1, dim))});

            const auto& corners = geometryHelper.getScvfCorners(scvfLocalIdx);
            scvfs_[scvfLocalIdx] = SubControlVolumeFace(
                corners,
                geometryHelper.normal(corners, localScvIndices),
                element,
                elementGeometry,
                scvfLocalIdx,
                std::move(localScvIndices),
                false
            );
        }

        // construct the sub control volume faces on the domain boundary
        for (const auto& intersection : intersections(gridGeometry().gridView(), element))
        {
            if (intersection.boundary() && !intersection.neighbor())
            {
                const auto isGeometry = intersection.geometry();
                hasBoundaryScvf_ = true;

                for (unsigned int isScvfLocalIdx = 0; isScvfLocalIdx < isGeometry.corners(); ++isScvfLocalIdx)
                {
                    // find the scv this scvf is connected to
                    const LocalIndexType insideScvIdx = static_cast<LocalIndexType>(refElement.subEntity(intersection.indexInInside(), 1, isScvfLocalIdx, dim));
                    std::vector<LocalIndexType> localScvIndices = {insideScvIdx, insideScvIdx};

                    scvfs_.emplace_back(
                        geometryHelper.getBoundaryScvfCorners(intersection.indexInInside(), isScvfLocalIdx),
                        intersection.centerUnitOuterNormal(),
                        intersection,
                        isGeometry,
                        isScvfLocalIdx,
                        scvfLocalIdx,
                        std::move(localScvIndices),
                        true
                    );

                    scvfBoundaryGeometryKeys_.emplace_back(std::array<LocalIndexType, 2>{{
                        static_cast<LocalIndexType>(intersection.indexInInside()),
                        static_cast<LocalIndexType>(isScvfLocalIdx)
                    }});

                    // increment local counter
                    scvfLocalIdx++;
                }
            }
        }
    }

    //! The bound element
    GridIndexType eIdx_;
    std::optional<Element> element_;

    //! The global geometry cache
    const GGCache* ggCache_;

    //! vectors to store the geometries locally after binding an element
    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;
    std::vector<std::array<LocalIndexType, 2>> scvfBoundaryGeometryKeys_;

    bool hasBoundaryScvf_ = false;
};

} // end namespace Dumux

#endif
