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
 * \brief Base class for the finite volume geometry vector for cell-centered TPFA models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element of the grid partition.
 */
#ifndef DUMUX_DISCRETIZATION_CCTPFA_GLOBAL_FVGEOMETRY_HH
#define DUMUX_DISCRETIZATION_CCTPFA_GLOBAL_FVGEOMETRY_HH

#include <dumux/common/elementmap.hh>
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/fvelementgeometry.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the finite volume geometry vector for cell-centered TPFA models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
template<class TypeTag, bool EnableGlobalFVGeometryCache>
class CCTpfaGlobalFVGeometry
{};

// specialization in case the FVElementGeometries are stored globally
template<class TypeTag>
class CCTpfaGlobalFVGeometry<TypeTag, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Element = typename GridView::template Codim<0>::Entity;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    //! The local class needs access to the scv, scvfs and the fv element geometry
    //! as they are globally cached
    friend typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

public:
    //! Constructor
    CCTpfaGlobalFVGeometry(const GridView& gridView)
    : gridView_(gridView), elementMap_(gridView) {}

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

    //! The total number of boundary sub control volume faces
    std::size_t numBoundaryScvf() const
    {
        return numBoundaryScvf_;
    }

    // Get an element from a sub control volume contained in it
    Element element(const SubControlVolume& scv) const
    { return elementMap_.element(scv.elementIndex()); }

    // Get an element from a global element index
    Element element(IndexType eIdx) const
    { return elementMap_.element(eIdx); }

    //! Return the gridView this global object lives on
    const GridView& gridView() const
    { return gridView_; }

    //! update all fvElementGeometries (do this again after grid adaption)
    void update(const Problem& problem)
    {
        problemPtr_ = &problem;

        // clear containers (necessary after grid refinement)
        scvs_.clear();
        scvfs_.clear();
        scvfIndicesOfScv_.clear();
        elementMap_.clear();

        // determine size of containers
        IndexType numScvs = gridView_.size(0);
        IndexType numScvf = 0;
        for (const auto& element : elements(gridView_))
            numScvf += element.subEntities(1);

        // reserve memory
        elementMap_.resize(numScvs);
        scvs_.resize(numScvs);
        scvfs_.reserve(numScvf);
        scvfIndicesOfScv_.resize(numScvs);

        // Build the scvs and scv faces
        IndexType scvfIdx = 0;
        numBoundaryScvf_ = 0;
        for (const auto& element : elements(gridView_))
        {
            auto eIdx = problem.elementMapper().index(element);
            scvs_[eIdx] = SubControlVolume(element.geometry(), eIdx);

            // fill the element map with seeds
            elementMap_[eIdx] = element.seed();

            // the element-wise index sets for finite volume geometry
            std::vector<IndexType> scvfsIndexSet;
            scvfsIndexSet.reserve(element.subEntities(1));

            // for network grids there might be multiple intersection with the same geometryInInside
            // we indentify those by the indexInInside for now (assumes conforming grids at branching facets)
            std::vector<std::vector<IndexType>> outsideIndices;
            if (dim < dimWorld)
            {
                outsideIndices.resize(element.subEntities(1));
                for (const auto& intersection : intersections(gridView_, element))
                {
                    if (intersection.neighbor())
                    {
                        const auto& outside = intersection.outside();
                        auto nIdx = problem.elementMapper().index(outside);
                        outsideIndices[intersection.indexInInside()].push_back(nIdx);
                    }

                }
            }

            for (const auto& intersection : intersections(gridView_, element))
            {
                // inner sub control volume faces
                if (intersection.neighbor())
                {
                    if (dim == dimWorld)
                    {
                        auto nIdx = problem.elementMapper().index(intersection.outside());
                        scvfs_.emplace_back(intersection,
                                            intersection.geometry(),
                                            scvfIdx,
                                            std::vector<IndexType>({eIdx, nIdx}));
                        scvfsIndexSet.push_back(scvfIdx++);
                    }
                    // this is for network grids
                    // (will be optimized away of dim == dimWorld)
                    else
                    {
                        auto indexInInside = intersection.indexInInside();
                        // check if we already handled this facet
                        if (outsideIndices[indexInInside].empty())
                            continue;
                        else
                        {
                            std::vector<IndexType> scvIndices({eIdx});
                            scvIndices.insert(scvIndices.end(), outsideIndices[indexInInside].begin(), outsideIndices[indexInInside].end());
                            scvfs_.emplace_back(intersection,
                                                intersection.geometry(),
                                                scvfIdx,
                                                scvIndices);
                            scvfsIndexSet.push_back(scvfIdx++);
                            outsideIndices[indexInInside].clear();
                        }
                    }
                }
                // boundary sub control volume faces
                else if (intersection.boundary())
                {
                    scvfs_.emplace_back(intersection,
                                        intersection.geometry(),
                                        scvfIdx,
                                        std::vector<IndexType>({eIdx, gridView_.size(0) + numBoundaryScvf_++})
                                        );
                    scvfsIndexSet.push_back(scvfIdx++);
                }
            }

            // Save the scvf indices belonging to this scv to build up fv element geometries fast
            scvfIndicesOfScv_[eIdx] = scvfsIndexSet;
        }

        // Make the flip index set for network and surface grids
        if (dim < dimWorld)
        {
            flipScvfIndices_.resize(scvfs_.size());
            for (auto&& scvf : scvfs_)
            {
                if (scvf.boundary())
                    continue;

                flipScvfIndices_[scvf.index()].resize(scvf.numOutsideScvs());
                const auto insideScvIdx = scvf.insideScvIdx();
                // check which outside scvf has the insideScvIdx index in its outsideScvIndices
                for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
                    flipScvfIndices_[scvf.index()][i] = findFlippedScvfIndex_(insideScvIdx, scvf.outsideScvIdx(i));
            }
        }
    }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline FVElementGeometry localView(const CCTpfaGlobalFVGeometry& global)
    { return FVElementGeometry(global); }

//private:

    //! Get a sub control volume with a global scv index
    const SubControlVolume& scv(IndexType scvIdx) const
    {
        return scvs_[scvIdx];
    }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& scvf(IndexType scvfIdx) const
    {
        return scvfs_[scvfIdx];
    }

    //! Get the scvf on the same face but from the other side
    //! Note that e.g. the normals might be different in the case of surface grids
    const SubControlVolumeFace& flipScvf(IndexType scvfIdx, unsigned int outsideScvfIdx = 0) const
    {
        return scvfs_[flipScvfIndices_[scvfIdx][outsideScvfIdx]];
    }

    //! Get the sub control volume face indices of an scv by global index
    const std::vector<IndexType>& scvfIndicesOfScv(IndexType scvIdx) const
    {
        return scvfIndicesOfScv_[scvIdx];
    }

private:
    // find the scvf that has insideScvIdx in its outsideScvIdx list and outsideScvIdx as its insideScvIdx
    IndexType findFlippedScvfIndex_(IndexType insideScvIdx, IndexType outsideScvIdx)
    {
        // go over all potential scvfs of the outside scv
        for (auto outsideScvfIndex : scvfIndicesOfScv_[outsideScvIdx])
        {
            const auto& outsideScvf = this->scvf(outsideScvfIndex);
            for (int j = 0; j < outsideScvf.numOutsideScvs(); ++j)
                if (outsideScvf.outsideScvIdx(j) == insideScvIdx)
                    return outsideScvf.index();
        }

        DUNE_THROW(Dune::InvalidStateException, "No flipped version of this scvf found!");
    }

    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;

    const GridView gridView_;
    Dumux::ElementMap<GridView> elementMap_;
    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;
    std::vector<std::vector<IndexType>> scvfIndicesOfScv_;
    IndexType numBoundaryScvf_;
    // needed for embedded surface and network grids (dim < dimWorld)
    std::vector<std::vector<IndexType>> flipScvfIndices_;
};

// specialization in case the FVElementGeometries are not stored
template<class TypeTag>
class CCTpfaGlobalFVGeometry<TypeTag, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Element = typename GridView::template Codim<0>::Entity;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    //! The local fvGeometry needs access to the problem
    friend typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

public:
    //! Constructor
    CCTpfaGlobalFVGeometry(const GridView& gridView)
    : gridView_(gridView), elementMap_(gridView) {}

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return numScvs_;
    }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    {
        return numScvf_;
    }

    //! The total number of boundary sub control volume faces
    std::size_t numBoundaryScvf() const
    {
        return numBoundaryScvf_;
    }

    // Get an element from a sub control volume contained in it
    Element element(const SubControlVolume& scv) const
    { return elementMap_.element(scv.elementIndex()); }

    // Get an element from a global element index
    Element element(IndexType eIdx) const
    { return elementMap_.element(eIdx); }

    //! Return the gridView this global object lives on
    const GridView& gridView() const
    { return gridView_; }

    //! update all fvElementGeometries (do this again after grid adaption)
    void update(const Problem& problem)
    {
        problemPtr_ = &problem;
        elementMap_.clear();

        // reserve memory or resize the containers
        numScvs_ = gridView_.size(0);
        numScvf_ = 0;
        numBoundaryScvf_ = 0;
        elementMap_.resize(numScvs_);
        scvfIndicesOfScv_.resize(numScvs_);
        neighborVolVarIndices_.resize(numScvs_);

        // Build the SCV and SCV face
        for (const auto& element : elements(gridView_))
        {
            auto eIdx = problem.elementMapper().index(element);

            // fill the element map with seeds
            elementMap_[eIdx] = element.seed();

            // the element-wise index sets for finite volume geometry
            auto numLocalFaces = element.subEntities(1);
            std::vector<IndexType> scvfsIndexSet;
            std::vector<std::vector<IndexType>> neighborVolVarIndexSet;
            scvfsIndexSet.reserve(numLocalFaces);
            neighborVolVarIndexSet.reserve(numLocalFaces);

            // for network grids there might be multiple intersection with the same geometryInInside
            // we indentify those by the indexInInside for now (assumes conforming grids at branching facets)
            std::vector<std::vector<IndexType>> outsideIndices;
            if (dim < dimWorld)
            {
                outsideIndices.resize(numLocalFaces);
                for (const auto& intersection : intersections(gridView_, element))
                {
                    if (intersection.neighbor())
                    {
                        const auto& outside = intersection.outside();
                        auto nIdx = problem.elementMapper().index(outside);
                        outsideIndices[intersection.indexInInside()].push_back(nIdx);
                    }

                }
            }

            for (const auto& intersection : intersections(gridView_, element))
            {
                // inner sub control volume faces
                if (intersection.neighbor())
                {
                    if (dim == dimWorld)
                    {
                        scvfsIndexSet.push_back(numScvf_++);
                        auto nIdx = problem.elementMapper().index(intersection.outside());
                        neighborVolVarIndexSet.push_back({nIdx});
                    }
                    // this is for network grids
                    // (will be optimized away of dim == dimWorld)
                    else
                    {
                        auto indexInInside = intersection.indexInInside();
                        // check if we already handled this facet
                        if (outsideIndices[indexInInside].empty())
                            continue;
                        else
                        {
                            scvfsIndexSet.push_back(numScvf_++);
                            neighborVolVarIndexSet.push_back(outsideIndices[indexInInside]);
                            outsideIndices[indexInInside].clear();
                        }
                    }
                }
                // boundary sub control volume faces
                else if (intersection.boundary())
                {
                    scvfsIndexSet.push_back(numScvf_++);
                    neighborVolVarIndexSet.push_back({numScvs_ + numBoundaryScvf_++});
                }
            }

            // store the sets of indices in the data container
            scvfIndicesOfScv_[eIdx] = scvfsIndexSet;
            neighborVolVarIndices_[eIdx] = neighborVolVarIndexSet;
        }
    }

    const std::vector<IndexType>& scvfIndicesOfScv(IndexType scvIdx) const
    { return scvfIndicesOfScv_[scvIdx]; }

    //! Return the neighbor volVar indices for all scvfs in the scv with index scvIdx
    const std::vector<std::vector<IndexType>>& neighborVolVarIndices(IndexType scvIdx) const
    { return neighborVolVarIndices_[scvIdx]; }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline FVElementGeometry localView(const CCTpfaGlobalFVGeometry& global)
    { return FVElementGeometry(global); }

private:

    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;

    const GridView gridView_;

    // Information on the global number of geometries
    IndexType numScvs_;
    IndexType numScvf_;
    IndexType numBoundaryScvf_;

    // vectors that store the global data
    Dumux::ElementMap<GridView> elementMap_;
    std::vector<std::vector<IndexType>> scvfIndicesOfScv_;
    std::vector<std::vector<std::vector<IndexType>>> neighborVolVarIndices_;
};

} // end namespace

#endif
