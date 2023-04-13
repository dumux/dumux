// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCTpfaDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for cell-centered TPFA models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element in the local scope we are restricting to, e.g. stencil or element.
 */
#ifndef DUMUX_DISCRETIZATION_CCTPFA_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_CCTPFA_FV_ELEMENT_GEOMETRY_HH

#include <optional>
#include <algorithm>
#include <array>
#include <vector>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dumux/common/indextraits.hh>
#include <dune/common/iteratorrange.hh>
#include <dumux/discretization/scvandscvfiterators.hh>

namespace Dumux {

namespace Detail::Tpfa {

template<class GridIndexType>
auto findLocalIndex(const GridIndexType idx,
                    const std::vector<GridIndexType>& indices)
{
    auto it = std::find(indices.begin(), indices.end(), idx);
    assert(it != indices.end() && "Could not find the scv/scvf! Make sure to properly bind this class!");
    return std::distance(indices.begin(), it);
}

} // end namespace Detail::Tpfa

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for cell-centered TPFA models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element in the local scope we are restricting to, e.g. stencil or element.
 * \tparam GG the finite volume grid geometry type
 * \tparam enableGridGeometryCache if the grid geometry is cached or not
 * \note This class is specialized for versions with and without caching the fv geometries on the grid view
 */
template<class GG, bool enableGridGeometryCache>
class CCTpfaFVElementGeometry;

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for cell-centered TPFA models
 *        Specialization for grid caching enabled
 * \note The finite volume geometries are stored in the corresponding GridGeometry
 */
template<class GG>
class CCTpfaFVElementGeometry<GG, true>
{
    using ThisType = CCTpfaFVElementGeometry<GG, true>;
    using GridView = typename GG::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;

public:
    //! export type of the element
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using GridGeometry = GG;

    //! the maximum number of scvs per element
    static constexpr std::size_t maxNumElementScvs = 1;
    //! the maximum number of scvfs per element (use cubes for maximum)
    static constexpr std::size_t maxNumElementScvfs = 2*GridView::dimension;

    //! Constructor
    CCTpfaFVElementGeometry(const GridGeometry& gridGeometry)
    : gridGeometryPtr_(&gridGeometry) {}

    //! Get an element sub control volume with a global scv index
    //! We separate element and neighbor scvs to speed up mapping
    const SubControlVolume& scv(GridIndexType scvIdx) const
    {
        return gridGeometry().scv(scvIdx);
    }

    //! Get an element sub control volume face with a global scvf index
    //! We separate element and neighbor scvfs to speed up mapping
    const SubControlVolumeFace& scvf(GridIndexType scvfIdx) const
    {
        return gridGeometry().scvf(scvfIdx);
    }

    //! Get the scvf on the same face but from the other side
    //! Note that e.g. the normals might be different in the case of surface grids
    const SubControlVolumeFace& flipScvf(GridIndexType scvfIdx, unsigned int outsideScvIdx = 0) const
    {
        return gridGeometry().flipScvf(scvfIdx, outsideScvIdx);
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element (not including neighbor scvs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange< ScvIterator<SubControlVolume, std::array<GridIndexType, 1>, ThisType> >
    scvs(const CCTpfaFVElementGeometry& fvGeometry)
    {
        using ScvIterator = Dumux::ScvIterator<SubControlVolume, std::array<GridIndexType, 1>, ThisType>;
        return Dune::IteratorRange<ScvIterator>(ScvIterator(fvGeometry.scvIndices_.begin(), fvGeometry),
                                                ScvIterator(fvGeometry.scvIndices_.end(), fvGeometry));
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element (not including neighbor scvfs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline Dune::IteratorRange< ScvfIterator<SubControlVolumeFace, std::vector<GridIndexType>, ThisType> >
    scvfs(const CCTpfaFVElementGeometry& fvGeometry)
    {
        const auto& g = fvGeometry.gridGeometry();
        const auto scvIdx = fvGeometry.scvIndices_[0];
        using ScvfIterator = Dumux::ScvfIterator<SubControlVolumeFace, std::vector<GridIndexType>, ThisType>;
        return Dune::IteratorRange<ScvfIterator>(ScvfIterator(g.scvfIndicesOfScv(scvIdx).begin(), fvGeometry),
                                                 ScvfIterator(g.scvfIndicesOfScv(scvIdx).end(), fvGeometry));
    }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScv() const
    {
        return scvIndices_.size();
    }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScvf() const
    {
        return gridGeometry().scvfIndicesOfScv(scvIndices_[0]).size();
    }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    CCTpfaFVElementGeometry bind(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    void bind(const Element& element) &
    {
        this->bindElement(element);
    }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bindElement(element);`
     */
    CCTpfaFVElementGeometry bindElement(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! Bind only element-local
    void bindElement(const Element& element) &
    {
        element_ = element;
        scvIndices_[0] = gridGeometry().elementMapper().index(*element_);
    }

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return static_cast<bool>(element_); }

    //! The bound element
    const Element& element() const
    { return *element_; }

    //! The global finite volume geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return *gridGeometryPtr_; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return gridGeometry().hasBoundaryScvf(scvIndices_[0]); }

    typename Element::Geometry geometry(const SubControlVolume& scv) const
    { return gridGeometryPtr_->element(scv.dofIndex()).geometry(); }

    typename GridView::Intersection::Geometry geometry(const SubControlVolumeFace& scvf) const
    {
        const auto element = gridGeometryPtr_->element(scvf.insideScvIdx());
        const auto& scvfIndices = gridGeometryPtr_->scvfIndicesOfScv(scvf.insideScvIdx());
        const LocalIndexType localScvfIdx = Detail::Tpfa::findLocalIndex(scvf.index(), scvfIndices);
        LocalIndexType localIdx = 0;
        for (const auto& intersection : intersections(gridGeometryPtr_->gridView(), element))
        {
            if (intersection.neighbor() || intersection.boundary())
            {
                if (localIdx == localScvfIdx)
                    return intersection.geometry();
                else
                    ++localIdx;
            }
        }

        DUNE_THROW(Dune::InvalidStateException, "Could not find scvf geometry");
    }

private:

    std::optional<Element> element_;
    std::array<GridIndexType, 1> scvIndices_;
    const GridGeometry* gridGeometryPtr_;
};

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for cell-centered TPFA models
 *        Specialization for grid caching disabled
 */
template<class GG>
class CCTpfaFVElementGeometry<GG, false>
{
    using ThisType = CCTpfaFVElementGeometry<GG, false>;
    using GridView = typename GG::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

public:
    //! export type of the element
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using GridGeometry = GG;
    //! the maximum number of scvs per element
    static constexpr std::size_t maxNumElementScvs = 1;
    //! the maximum number of scvfs per element (use cubes for maximum)
    static constexpr std::size_t maxNumElementScvfs = 2*dim;

    //! Constructor
    CCTpfaFVElementGeometry(const GridGeometry& gridGeometry)
    : gridGeometryPtr_(&gridGeometry) {}

    //! Get an element sub control volume with a global scv index
    //! We separate element and neighbor scvs to speed up mapping
    const SubControlVolume& scv(GridIndexType scvIdx) const
    {
        if (scvIdx == scvIndices_[0])
            return scvs_[0];
        else
            return neighborScvs_[Detail::Tpfa::findLocalIndex(scvIdx, neighborScvIndices_)];
    }

    //! Get an element sub control volume face with a global scvf index
    //! We separate element and neighbor scvfs to speed up mapping
    const SubControlVolumeFace& scvf(GridIndexType scvfIdx) const
    {
        auto it = std::find(scvfIndices_.begin(), scvfIndices_.end(), scvfIdx);
        if (it != scvfIndices_.end())
            return scvfs_[std::distance(scvfIndices_.begin(), it)];
        else
            return neighborScvfs_[Detail::Tpfa::findLocalIndex(scvfIdx, neighborScvfIndices_)];
    }

    //! Get the scvf on the same face but from the other side
    //! Note that e.g. the normals might be different in the case of surface grids
    const SubControlVolumeFace& flipScvf(GridIndexType scvfIdx, unsigned int outsideScvIdx = 0) const
    {
        auto it = std::find(scvfIndices_.begin(), scvfIndices_.end(), scvfIdx);
        if (it != scvfIndices_.end())
        {
            const auto localScvfIdx = std::distance(scvfIndices_.begin(), it);
            return neighborScvfs_[flippedScvfIndices_[localScvfIdx][outsideScvIdx]];
        }
        else
        {
            const auto localScvfIdx = Detail::Tpfa::findLocalIndex(scvfIdx, neighborScvfIndices_);
            const auto localFlippedIndex = flippedNeighborScvfIndices_[localScvfIdx][outsideScvIdx];
            if (localFlippedIndex < scvfs_.size())
                return scvfs_[localFlippedIndex];
            else
                return neighborScvfs_[localFlippedIndex - scvfs_.size()];
        }
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element (not including neighbor scvs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::array<SubControlVolume, 1>::const_iterator>
    scvs(const ThisType& g)
    {
        using IteratorType = typename std::array<SubControlVolume, 1>::const_iterator;
        return Dune::IteratorRange<IteratorType>(g.scvs_.begin(), g.scvs_.end());
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element (not including neighbor scvfs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolumeFace>::const_iterator>
    scvfs(const ThisType& g)
    {
        using IteratorType = typename std::vector<SubControlVolumeFace>::const_iterator;
        return Dune::IteratorRange<IteratorType>(g.scvfs_.begin(), g.scvfs_.end());
    }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScv() const
    { return scvs_.size(); }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScvf() const
    { return scvfs_.size(); }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    CCTpfaFVElementGeometry bind(const Element& element) &&
    {
        this->bind_(element);
        return std::move(*this);
    }

    void bind(const Element& element) &
    {
        this->bind_(element);
    }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bindElement(element);`
     */
    CCTpfaFVElementGeometry bindElement(const Element& element) &&
    {
        this->bindElement_(element);
        return std::move(*this);
    }

    void bindElement(const Element& element) &
    {
        this->bindElement_(element);
    }

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return static_cast<bool>(element_); }

    //! The bound element
    const Element& element() const
    { return *element_; }

    //! The global finite volume geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return *gridGeometryPtr_; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return hasBoundaryScvf_; }

    typename Element::Geometry geometry(const SubControlVolume& scv) const
    { return gridGeometryPtr_->element(scv.dofIndex()).geometry(); }

    typename GridView::Intersection::Geometry geometry(const SubControlVolumeFace& scvf) const
    {
        const auto element = gridGeometryPtr_->element(scvf.insideScvIdx());
        const auto& scvfIndices = gridGeometryPtr_->scvfIndicesOfScv(scvf.insideScvIdx());
        const LocalIndexType localScvfIdx = Detail::Tpfa::findLocalIndex(scvf.index(), scvfIndices);
        LocalIndexType localIdx = 0;
        for (const auto& intersection : intersections(gridGeometryPtr_->gridView(), element))
        {
            if (intersection.neighbor() || intersection.boundary())
            {
                if (localIdx == localScvfIdx)
                    return intersection.geometry();
                else
                    ++localIdx;
            }
        }

        DUNE_THROW(Dune::InvalidStateException, "Could not find scvf geometry");
    }

private:
    //! Binding of an element preparing the geometries of the whole stencil
    //! called by the local jacobian to prepare element assembly
    void bind_(const Element& element)
    {
        bindElement_(element);

        neighborScvs_.reserve(element.subEntities(1));
        neighborScvfIndices_.reserve(element.subEntities(1));
        neighborScvfs_.reserve(element.subEntities(1));

        std::vector<GridIndexType> handledNeighbors;
        handledNeighbors.reserve(element.subEntities(1));

        for (const auto& intersection : intersections(gridGeometry().gridView(), element))
        {
            // for inner intersections and periodic (according to grid interface) intersections make neighbor geometry
            if (intersection.neighbor())
            {
                const auto outside = intersection.outside();
                const auto outsideIdx = gridGeometry().elementMapper().index(outside);

                // make outside geometries only if not done yet (could happen on non-conforming grids)
                if ( std::find(handledNeighbors.begin(), handledNeighbors.end(), outsideIdx) == handledNeighbors.end() )
                {
                    makeNeighborGeometries(outside, outsideIdx);
                    handledNeighbors.push_back(outsideIdx);
                }
            }
        }

        // build flip index set for network, surface, and periodic grids
        if (dim < dimWorld || gridGeometry().isPeriodic())
        {
            flippedScvfIndices_.resize(scvfs_.size());
            for (unsigned int localScvfIdx = 0; localScvfIdx < scvfs_.size(); ++localScvfIdx)
            {
                const auto& scvf = scvfs_[localScvfIdx];
                if (scvf.boundary())
                    continue;

                flippedScvfIndices_[localScvfIdx].resize(scvf.numOutsideScvs());
                for (unsigned int localOutsideScvIdx = 0; localOutsideScvIdx < scvf.numOutsideScvs(); ++localOutsideScvIdx)
                {
                    const auto globalOutsideScvIdx = scvf.outsideScvIdx(localOutsideScvIdx);
                    for (unsigned int localNeighborScvfIdx = 0; localNeighborScvfIdx < neighborScvfs_.size(); ++localNeighborScvfIdx)
                    {
                        if (neighborScvfs_[localNeighborScvfIdx].insideScvIdx() == globalOutsideScvIdx)
                        {
                            flippedScvfIndices_[localScvfIdx][localOutsideScvIdx] = localNeighborScvfIdx;
                            break;
                        }
                    }
                }
            }

            flippedNeighborScvfIndices_.resize(neighborScvfs_.size());
            for (unsigned int localScvfIdx = 0; localScvfIdx < neighborScvfs_.size(); ++localScvfIdx)
            {
                const auto& neighborScvf = neighborScvfs_[localScvfIdx];
                flippedNeighborScvfIndices_[localScvfIdx].resize(neighborScvf.numOutsideScvs());
                for (unsigned int localOutsideScvIdx = 0; localOutsideScvIdx < neighborScvf.numOutsideScvs(); ++localOutsideScvIdx)
                {
                    flippedNeighborScvfIndices_[localScvfIdx][localOutsideScvIdx] = findFlippedScvfIndex_(neighborScvf.insideScvIdx(), neighborScvf.outsideScvIdx(localOutsideScvIdx));
                }
            }
        }
    }

    //! Binding of an element preparing the geometries only inside the element
    void bindElement_(const Element& element)
    {
        clear();
        element_ = element;
        scvfs_.reserve(element.subEntities(1));
        scvfIndices_.reserve(element.subEntities(1));
        makeElementGeometries(element);
    }

    GridIndexType findFlippedScvfIndex_(GridIndexType insideScvIdx, GridIndexType globalOutsideScvIdx)
    {
        for (unsigned int localNeighborScvfIdx = 0; localNeighborScvfIdx < neighborScvfs_.size(); ++localNeighborScvfIdx)
        {
            if (neighborScvfs_[localNeighborScvfIdx].insideScvIdx() == globalOutsideScvIdx)
            {
                return scvfs_.size() + localNeighborScvfIdx;
            }
        }

        // go over all potential scvfs of the outside scv
        for (unsigned int localOutsideScvfIdx = 0; localOutsideScvfIdx < scvfs_.size(); ++localOutsideScvfIdx)
        {
            const auto& outsideScvf = scvfs_[localOutsideScvfIdx];
            for (unsigned int j = 0; j < outsideScvf.numOutsideScvs(); ++j)
            {
                if (outsideScvf.outsideScvIdx(j) == insideScvIdx)
                {
                    return localOutsideScvfIdx;
                }
            }
        }

        DUNE_THROW(Dune::InvalidStateException, "No flipped version of this scvf found!");
    }

    //! create scvs and scvfs of the bound element
    void makeElementGeometries(const Element& element)
    {
        using ScvfGridIndexStorage = typename SubControlVolumeFace::Traits::GridIndexStorage;

        const auto eIdx = gridGeometry().elementMapper().index(element);
        scvs_[0] = SubControlVolume(element.geometry(), eIdx);
        scvIndices_[0] = eIdx;

        const auto& scvFaceIndices = gridGeometry().scvfIndicesOfScv(eIdx);
        const auto& neighborVolVarIndices = gridGeometry().neighborVolVarIndices(eIdx);

        // for network grids there might be multiple intersection with the same geometryInInside
        // we identify those by the indexInInside for now (assumes conforming grids at branching facets)
        // here we keep track of them
        std::vector<bool> handledScvf;
        if (dim < dimWorld)
            handledScvf.resize(element.subEntities(1), false);

        int scvfCounter = 0;
        for (const auto& intersection : intersections(gridGeometry().gridView(), element))
        {
            if (dim < dimWorld)
                if (handledScvf[intersection.indexInInside()])
                    continue;

            const auto& scvfNeighborVolVarIndices = neighborVolVarIndices[scvfCounter];
            if (intersection.neighbor() || intersection.boundary())
            {
                ScvfGridIndexStorage scvIndices;
                scvIndices.resize(scvfNeighborVolVarIndices.size() + 1);
                scvIndices[0] = eIdx;
                std::copy(scvfNeighborVolVarIndices.begin(), scvfNeighborVolVarIndices.end(), scvIndices.begin()+1);

                const bool onBoundary = intersection.boundary() && !intersection.neighbor();
                hasBoundaryScvf_ = (hasBoundaryScvf_ || onBoundary);

                scvfs_.emplace_back(intersection,
                                    intersection.geometry(),
                                    scvFaceIndices[scvfCounter],
                                    scvIndices,
                                    onBoundary);
                scvfIndices_.emplace_back(scvFaceIndices[scvfCounter]);
                scvfCounter++;

                // for surface and network grids mark that we handled this face
                if (dim < dimWorld)
                    handledScvf[intersection.indexInInside()] = true;
            }
        }
    }

    //! create the necessary scvs and scvfs of the neighbor elements to the bound elements
    void makeNeighborGeometries(const Element& element, const GridIndexType eIdx)
    {
        using ScvfGridIndexStorage = typename SubControlVolumeFace::Traits::GridIndexStorage;

        // create the neighbor scv
        neighborScvs_.emplace_back(element.geometry(), eIdx);
        neighborScvIndices_.push_back(eIdx);

        const auto& scvFaceIndices = gridGeometry().scvfIndicesOfScv(eIdx);
        const auto& neighborVolVarIndices = gridGeometry().neighborVolVarIndices(eIdx);

        // for network grids there might be multiple intersection with the same geometryInInside
        // we identify those by the indexInInside for now (assumes conforming grids at branching facets)
        // here we keep track of them
        std::vector<bool> handledScvf;
        if (dim < dimWorld)
            handledScvf.resize(element.subEntities(1), false);

        int scvfCounter = 0;
        for (const auto& intersection : intersections(gridGeometry().gridView(), element))
        {
            if (dim < dimWorld)
                if (handledScvf[intersection.indexInInside()])
                    continue;

            if (intersection.neighbor())
            {
                // this catches inner and periodic scvfs
                const auto& scvfNeighborVolVarIndices = neighborVolVarIndices[scvfCounter];
                if (scvfNeighborVolVarIndices[0] < gridGeometry().gridView().size(0))
                {
                    // only create subcontrol faces where the outside element is the bound element
                    if (dim == dimWorld)
                    {
                        if (scvfNeighborVolVarIndices[0] == gridGeometry().elementMapper().index(*element_))
                        {
                            ScvfGridIndexStorage scvIndices({eIdx, scvfNeighborVolVarIndices[0]});
                            neighborScvfs_.emplace_back(intersection,
                                                        intersection.geometry(),
                                                        scvFaceIndices[scvfCounter],
                                                        scvIndices,
                                                        false);

                            neighborScvfIndices_.push_back(scvFaceIndices[scvfCounter]);
                        }
                    }
                    // for network grids we can't use the intersection.outside() index as we can't assure that the
                    // first intersection with this indexInInside is the one that has our bound element as outside
                    // instead we check if the bound element's index is in the outsideScvIndices of the candidate scvf
                    // (will be optimized away for dim == dimWorld)
                    else
                    {
                        for (unsigned outsideScvIdx = 0; outsideScvIdx < scvfNeighborVolVarIndices.size(); ++outsideScvIdx)
                        {
                            if (scvfNeighborVolVarIndices[outsideScvIdx] == gridGeometry().elementMapper().index(*element_))
                            {
                                ScvfGridIndexStorage scvIndices;
                                scvIndices.resize(scvfNeighborVolVarIndices.size() + 1);
                                scvIndices[0] = eIdx;
                                std::copy(scvfNeighborVolVarIndices.begin(), scvfNeighborVolVarIndices.end(), scvIndices.begin()+1);
                                neighborScvfs_.emplace_back(intersection,
                                                            intersection.geometry(),
                                                            scvFaceIndices[scvfCounter],
                                                            scvIndices,
                                                            false);

                                neighborScvfIndices_.push_back(scvFaceIndices[scvfCounter]);
                                break;
                            }
                        }
                    }

                    // for surface and network grids mark that we handled this face
                    if (dim < dimWorld)
                        handledScvf[intersection.indexInInside()] = true;
                    scvfCounter++;
                }
            }

            // only increase counter for boundary intersections
            // (exclude periodic boundaries according to dune grid interface, they have been handled in neighbor==true)
            else if (intersection.boundary())
            {
                if (dim < dimWorld)
                    handledScvf[intersection.indexInInside()] = true;
                scvfCounter++;
            }
        }
    }

    //! Clear all local data
    void clear()
    {
        scvfIndices_.clear();
        scvfs_.clear();
        flippedScvfIndices_.clear();

        neighborScvIndices_.clear();
        neighborScvfIndices_.clear();
        neighborScvs_.clear();
        neighborScvfs_.clear();
        flippedNeighborScvfIndices_.clear();

        hasBoundaryScvf_ = false;
    }

    std::optional<Element> element_; //!< the element to which this fvgeometry is bound
    const GridGeometry* gridGeometryPtr_;  //!< the grid fvgeometry

    // local storage after binding an element
    std::array<GridIndexType, 1> scvIndices_;
    std::array<SubControlVolume, 1> scvs_;

    std::vector<GridIndexType> scvfIndices_;
    std::vector<SubControlVolumeFace> scvfs_;
    std::vector<std::vector<GridIndexType>> flippedScvfIndices_;

    std::vector<GridIndexType> neighborScvIndices_;
    std::vector<SubControlVolume> neighborScvs_;

    std::vector<GridIndexType> neighborScvfIndices_;
    std::vector<SubControlVolumeFace> neighborScvfs_;
    std::vector<std::vector<GridIndexType>> flippedNeighborScvfIndices_;

    bool hasBoundaryScvf_ = false;
};

} // end namespace Dumux

#endif
