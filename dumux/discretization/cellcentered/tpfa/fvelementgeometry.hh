// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup CCTpfaDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for cell-centered TPFA models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element in the local scope we are restricting to, e.g. stencil or element.
 */
#ifndef DUMUX_DISCRETIZATION_CCTPFA_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_CCTPFA_FV_ELEMENT_GEOMETRY_HH

#include <dune/common/exceptions.hh>
#include <dune/common/iteratorrange.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/scvandscvfiterators.hh>

namespace Dumux
{

//! forward declaration of the global finite volume geometry
template<class TypeTag, bool EnableFVGridGeometryCache>
class CCTpfaFVGridGeometry;

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for cell-centered TPFA models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element in the local scope we are restricting to, e.g. stencil or element.
 * \note This class is specialized for versions with and without caching the fv geometries on the grid view
 */
template<class TypeTag, bool EnableFVGridGeometryCache>
class CCTpfaFVElementGeometry
{};

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for cell-centered TPFA models
 *        Specialization for grid caching enabled
 * \note The finite volume geometries are stored in the corresponding FVGridGeometry
 */
template<class TypeTag>
class CCTpfaFVElementGeometry<TypeTag, true>
{
    using ThisType = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! export type of subcontrol volume
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    //! export type of finite volume grid geometry
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    //! the maximum number of scvs per element
    static constexpr std::size_t maxNumElementScvs = 1;
    //! the maximum number of scvfs per element (use cubes for maximum)
    static constexpr std::size_t maxNumElementScvfs = 2*GridView::dimension;

    //! Constructor
    CCTpfaFVElementGeometry(const FVGridGeometry& fvGridGeometry)
    : fvGridGeometryPtr_(&fvGridGeometry) {}

    //! Get an elment sub control volume with a global scv index
    //! We separate element and neighbor scvs to speed up mapping
    const SubControlVolume& scv(IndexType scvIdx) const
    {
        return fvGridGeometry().scv(scvIdx);
    }

    //! Get an element sub control volume face with a global scvf index
    //! We separate element and neighbor scvfs to speed up mapping
    const SubControlVolumeFace& scvf(IndexType scvfIdx) const
    {
        return fvGridGeometry().scvf(scvfIdx);
    }

    //! Get the scvf on the same face but from the other side
    //! Note that e.g. the normals might be different in the case of surface grids
    const SubControlVolumeFace& flipScvf(IndexType scvfIdx, unsigned int outsideScvIdx = 0) const
    {
        return fvGridGeometry().flipScvf(scvfIdx, outsideScvIdx);
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element (not including neighbor scvs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange< ScvIterator<SubControlVolume, std::array<IndexType, 1>, ThisType> >
    scvs(const CCTpfaFVElementGeometry& fvGeometry)
    {
        using ScvIterator = Dumux::ScvIterator<SubControlVolume, std::array<IndexType, 1>, ThisType>;
        return Dune::IteratorRange<ScvIterator>(ScvIterator(fvGeometry.scvIndices_.begin(), fvGeometry),
                                                ScvIterator(fvGeometry.scvIndices_.end(), fvGeometry));
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element (not including neighbor scvfs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline Dune::IteratorRange< ScvfIterator<SubControlVolumeFace, std::vector<IndexType>, ThisType> >
    scvfs(const CCTpfaFVElementGeometry& fvGeometry)
    {
        const auto& g = fvGeometry.fvGridGeometry();
        const auto scvIdx = fvGeometry.scvIndices_[0];
        using ScvfIterator = Dumux::ScvfIterator<SubControlVolumeFace, std::vector<IndexType>, ThisType>;
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
        return fvGridGeometry().scvfIndicesOfScv(scvIndices_[0]).size();
    }

    //! Binding of an element, called by the local jacobian to prepare element assembly
    void bind(const Element& element)
    {
        this->bindElement(element);
    }

    //! Bind only element-local
    void bindElement(const Element& element)
    {
        elementPtr_ = &element;
        scvIndices_[0] = fvGridGeometry().elementMapper().index(*elementPtr_);
    }

    //! The global finite volume geometry we are a restriction of
    const FVGridGeometry& fvGridGeometry() const
    { return *fvGridGeometryPtr_; }

private:

    const Element* elementPtr_;
    std::array<IndexType, 1> scvIndices_;
    const FVGridGeometry* fvGridGeometryPtr_;
};

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for cell-centered TPFA models
 *        Specialization for grid caching disabled
 */
template<class TypeTag>
class CCTpfaFVElementGeometry<TypeTag, false>
{
    using ThisType = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

public:
    //! export type of subcontrol volume
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    //! export type of finite volume grid geometry
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    //! the maximum number of scvs per element
    static constexpr std::size_t maxNumElementScvs = 1;
    //! the maximum number of scvfs per element (use cubes for maximum)
    static constexpr std::size_t maxNumElementScvfs = 2*dim;

    //! Constructor
    CCTpfaFVElementGeometry(const FVGridGeometry& fvGridGeometry)
    : fvGridGeometryPtr_(&fvGridGeometry) {}

    //! Get an elment sub control volume with a global scv index
    //! We separate element and neighbor scvs to speed up mapping
    const SubControlVolume& scv(IndexType scvIdx) const
    {
        if (scvIdx == scvIndices_[0])
            return scvs_[0];
        else
            return neighborScvs_[findLocalIndex(scvIdx, neighborScvIndices_)];
    }

    //! Get an element sub control volume face with a global scvf index
    //! We separate element and neighbor scvfs to speed up mapping
    const SubControlVolumeFace& scvf(IndexType scvfIdx) const
    {
        auto it = std::find(scvfIndices_.begin(), scvfIndices_.end(), scvfIdx);
        if (it != scvfIndices_.end())
            return scvfs_[std::distance(scvfIndices_.begin(), it)];
        else
            return neighborScvfs_[findLocalIndex(scvfIdx, neighborScvfIndices_)];
    }

    //! Get the scvf on the same face but from the other side
    //! Note that e.g. the normals might be different in the case of surface grids
    const SubControlVolumeFace& flipScvf(IndexType scvfIdx, unsigned int outsideScvIdx = 0) const
    {
        auto it = std::find(scvfIndices_.begin(), scvfIndices_.end(), scvfIdx);
        if (it != scvfIndices_.end())
        {
            const auto localScvfIdx = std::distance(scvfIndices_.begin(), it);
            return neighborScvfs_[flippedScvfIndices_[localScvfIdx][outsideScvIdx]];
        }
        else
        {
            const auto localScvfIdx = findLocalIndex(scvfIdx, neighborScvfIndices_);
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

    //! Binding of an element preparing the geometries of the whole stencil
    //! called by the local jacobian to prepare element assembly
    void bind(const Element& element)
    {
        bindElement(element);

        neighborScvs_.reserve(element.subEntities(1));
        neighborScvfIndices_.reserve(element.subEntities(1));
        neighborScvfs_.reserve(element.subEntities(1));

        std::vector<IndexType> handledNeighbors;
        handledNeighbors.reserve(element.subEntities(1));
        for (const auto& intersection : intersections(fvGridGeometry().gridView(), element))
        {
            if (intersection.neighbor())
            {
                const auto outside = intersection.outside();
                const auto outsideIdx = fvGridGeometry().elementMapper().index(outside);

                // make outside geometries only if not done yet (could happen on non-conforming grids)
                if ( std::find(handledNeighbors.begin(), handledNeighbors.end(), outsideIdx) == handledNeighbors.end() )
                {
                    makeNeighborGeometries(outside, outsideIdx);
                    handledNeighbors.push_back(outsideIdx);
                }
            }
        }

        if (dim < dimWorld)
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

        //! TODO Check if user added additional DOF dependencies, i.e. the residual of DOF globalI depends
        //! on additional DOFs not included in the discretization schemes' occupation pattern
        // const auto globalI = fvGridGeometry().elementMapper().index(element);
        // const auto& additionalDofDependencies = problem.getAdditionalDofDependencies(globalI);
        // if (!additionalDofDependencies.empty())
        // {
        //     neighborScvs_.reserve(neighborScvs_.size() + additionalDofDependencies.size());
        //     neighborScvIndices_.reserve(neighborScvIndices_.size() + additionalDofDependencies.size());
        //     for (auto globalJ : additionalDofDependencies)
        //     {
        //         neighborScvs_.emplace_back(fvGridGeometry().element(globalJ).geometry(), globalJ);
        //         neighborScvIndices_.emplace_back(globalJ);
        //     }
        // }
    }

    //! Binding of an element preparing the geometries only inside the element
    void bindElement(const Element& element)
    {
        clear();
        elementPtr_ = &element;
        scvfs_.reserve(element.subEntities(1));
        scvfIndices_.reserve(element.subEntities(1));
        makeElementGeometries(element);
    }

    //! The global finite volume geometry we are a restriction of
    const FVGridGeometry& fvGridGeometry() const
    { return *fvGridGeometryPtr_; }

private:

    IndexType findFlippedScvfIndex_(IndexType insideScvIdx, IndexType globalOutsideScvIdx)
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
        const auto eIdx = fvGridGeometry().elementMapper().index(element);
        scvs_[0] = SubControlVolume(element.geometry(), eIdx);
        scvIndices_[0] = eIdx;

        const auto& scvFaceIndices = fvGridGeometry().scvfIndicesOfScv(eIdx);
        const auto& neighborVolVarIndices = fvGridGeometry().neighborVolVarIndices(eIdx);

        // for network grids there might be multiple intersection with the same geometryInInside
        // we indentify those by the indexInInside for now (assumes conforming grids at branching facets)
        // here we keep track of them
        std::vector<bool> handledScvf;
        if (dim < dimWorld)
            handledScvf.resize(element.subEntities(1), false);

        int scvfCounter = 0;
        for (const auto& intersection : intersections(fvGridGeometry().gridView(), element))
        {
            if (dim < dimWorld)
                if (handledScvf[intersection.indexInInside()])
                    continue;

            if (intersection.neighbor() || intersection.boundary())
            {
                std::vector<IndexType> scvIndices({eIdx});
                scvIndices.insert(scvIndices.end(), neighborVolVarIndices[scvfCounter].begin(), neighborVolVarIndices[scvfCounter].end());
                scvfs_.emplace_back(intersection,
                                    intersection.geometry(),
                                    scvFaceIndices[scvfCounter],
                                    scvIndices,
                                    intersection.boundary());
                scvfIndices_.emplace_back(scvFaceIndices[scvfCounter]);
                scvfCounter++;

                // for surface and network grids mark that we handled this face
                if (dim < dimWorld)
                    handledScvf[intersection.indexInInside()] = true;
            }
        }
    }

    //! create the necessary scvs and scvfs of the neighbor elements to the bound elements
    void makeNeighborGeometries(const Element& element, const IndexType eIdx)
    {
        // create the neighbor scv
        neighborScvs_.emplace_back(element.geometry(), eIdx);
        neighborScvIndices_.push_back(eIdx);

        const auto& scvFaceIndices = fvGridGeometry().scvfIndicesOfScv(eIdx);
        const auto& neighborVolVarIndices = fvGridGeometry().neighborVolVarIndices(eIdx);

        // for network grids there might be multiple intersection with the same geometryInInside
        // we indentify those by the indexInInside for now (assumes conforming grids at branching facets)
        // here we keep track of them
        std::vector<bool> handledScvf;
        if (dim < dimWorld)
            handledScvf.resize(element.subEntities(1), false);

        int scvfCounter = 0;
        for (const auto& intersection : intersections(fvGridGeometry().gridView(), element))
        {
            if (dim < dimWorld)
                if (handledScvf[intersection.indexInInside()])
                    continue;

            if (intersection.neighbor())
            {
                // only create subcontrol faces where the outside element is the bound element
                if (dim == dimWorld)
                {
                    if (intersection.outside() == *elementPtr_)
                    {
                        std::vector<IndexType> scvIndices({eIdx});
                        scvIndices.insert(scvIndices.end(), neighborVolVarIndices[scvfCounter].begin(), neighborVolVarIndices[scvfCounter].end());
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
                    for (unsigned outsideScvIdx = 0; outsideScvIdx < neighborVolVarIndices[scvfCounter].size(); ++outsideScvIdx)
                    {
                        if (neighborVolVarIndices[scvfCounter][outsideScvIdx] == fvGridGeometry().elementMapper().index(*elementPtr_))
                        {
                            std::vector<IndexType> scvIndices({eIdx});
                            scvIndices.insert(scvIndices.end(), neighborVolVarIndices[scvfCounter].begin(), neighborVolVarIndices[scvfCounter].end());
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
            else if (intersection.boundary())
            {
                // for surface and network grids mark that we handled this face
                if (dim < dimWorld)
                    handledScvf[intersection.indexInInside()] = true;
                scvfCounter++;
            }
        }
    }

    const IndexType findLocalIndex(const IndexType idx,
                                   const std::vector<IndexType>& indices) const
    {
        auto it = std::find(indices.begin(), indices.end(), idx);
        assert(it != indices.end() && "Could not find the scv/scvf! Make sure to properly bind this class!");
        return std::distance(indices.begin(), it);

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
    }

    const Element* elementPtr_; //!< the element to which this fvgeometry is bound
    const FVGridGeometry* fvGridGeometryPtr_;  //!< the grid fvgeometry

    // local storage after binding an element
    std::array<IndexType, 1> scvIndices_;
    std::array<SubControlVolume, 1> scvs_;

    std::vector<IndexType> scvfIndices_;
    std::vector<SubControlVolumeFace> scvfs_;
    std::vector<std::vector<IndexType>> flippedScvfIndices_;

    std::vector<IndexType> neighborScvIndices_;
    std::vector<SubControlVolume> neighborScvs_;

    std::vector<IndexType> neighborScvfIndices_;
    std::vector<SubControlVolumeFace> neighborScvfs_;
    std::vector<std::vector<IndexType>> flippedNeighborScvfIndices_;
};

} // end namespace Dumux

#endif
