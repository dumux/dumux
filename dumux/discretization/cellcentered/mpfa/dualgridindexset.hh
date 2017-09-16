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
 * \brief Base class for the index sets of the dual grid in mpfa schemes.
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_DUALGRID_INDEXSET_BASE_HH
#define DUMUX_DISCRETIZATION_MPFA_DUALGRID_INDEXSET_BASE_HH

#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/fvgridgeometry.hh>

namespace Dumux
{

/*!
 * \ingroup Mpfa
 * \brief Nodal index set for the dual grid of mpfa schemes.
 */
template<class TypeTag>
class DualGridNodalIndexSet
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using PrimaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
    using LocalIndexType = typename PrimaryInteractionVolume::Traits::LocalIndexType;
    using LocalIndexContainer = typename PrimaryInteractionVolume::Traits::DynamicLocalIndexContainer;
    using GlobalIndexContainer = typename PrimaryInteractionVolume::Traits::DynamicGlobalIndexContainer;

    static const int dim = GridView::dimension;

public:
    //! constructor
    DualGridNodalIndexSet() : numBoundaryScvfs_(0) {}

    // Inserts data for a given scvf. This should only called once per scvf!
    void insert(const SubControlVolumeFace& scvf)
    {
        insert(scvf.boundary(),
               scvf.index(),
               scvf.insideScvIdx(),
               scvf.outsideScvIndices());
    }

    // Inserts data for a given scvf. This should only called once per scvf!
    void insert(const bool boundary,
                const IndexType scvfIdx,
                const IndexType insideScvIdx,
                const std::vector<IndexType>& outsideScvIndices)
    {
        // check if it this really hasn't been called yet
        assert( std::find_if( scvfIndices_.begin(),
                              scvfIndices_.end(),
                              scvfIdx
                            ) == scvfIndices_.end() && "scvf has already been inserted!");

        // the local index of the scvf data about to be inserted
        const LocalIndexType curScvfLocalIdx = scvfIndices_.size();

        // add global scvf data
        GlobalIndexContainer scvIndices;
        scvIndices.reserve(outsideScvIndices.size()+1);
        scvIndices.push_back(insideScvIdx);
        scvIndices.insert(scvIndices.end(), outsideScvIndices.begin(), outsideScvIndices.end());

        // if scvf is on boundary, increase counter
        if (boundary)
            numBoundaryScvfs_++;

        scvfIndices_.push_back(scvfIdx);
        scvfIsOnBoundary_.push_back(boundary);
        scvfNeighborScvIndicesGlobal_.emplace_back( std::move(scvIndices) );

        // if entry for the inside scv exists append scvf local index, create otherwise
        auto it = std::find( scvIndices_.begin(), scvIndices_.end(), insideScvIdx );
        if (it != scvIndices_.end())
        {
            const auto localScvIdx = std::distance(scvIndices_.begin(), it);
            localScvfIndicesInScv_[localScvIdx].push_back(curScvfLocalIdx);
        }
        else
        {
            LocalIndexContainer localScvfIndices;
            localScvfIndices.reserve(dim);
            localScvfIndices.push_back(curScvfLocalIdx);
            localScvfIndicesInScv_.emplace_back( std::move(localScvfIndices) );
            scvIndices_.push_back(insideScvIdx);
        }
    }

    // This should be called AFTER (!) all scvf data has been inserted
    void makeLocal()
    {
        // make sure the created index set makes sense so far
        assert(checkGlobalIndexSet_());

        // compute local neighboring scv indices for the scvfs
        scvfNeighborScvIndices_.resize(numScvfs());
        for (unsigned int i = 0; i < numScvfs(); ++i)
        {
            const auto& neighborsGlobal = scvfNeighborScvIndicesGlobal_[i];
            const auto numNeighbors = scvfIsOnBoundary_[i] ? 1 : neighborsGlobal.size();

            scvfNeighborScvIndices_[i].resize(numNeighbors);
            for (unsigned int j = 0; j < numNeighbors; ++j)
                scvfNeighborScvIndices_[i][j] = findLocalScvIdx_(neighborsGlobal[j]);
        }

        // delete global neighboring scv data (not used anymore)
        scvfNeighborScvIndicesGlobal_.clear();
    }

    //! returns the number of scvs
    std::size_t numScvs() const
    { return scvIndices_.size(); }

    //! returns the number of scvfs
    std::size_t numScvfs() const
    { return scvfIndices_.size(); }

    //! returns the number of boundary scvfs
    std::size_t numBoundaryScvfs() const
    { return numBoundaryScvfs_; }

    //! returns the global scv indices connected to this dual grid node
    const GlobalIndexContainer& globalScvIndices() const
    { return scvIndices_; }

    //! returns the global scvf indices connected to this dual grid node
    const GlobalIndexContainer& globalScvfIndices() const
    { return scvfIndices_; }

    //! returns a global scv idx for a given local scv index
    IndexType scvIdxGlobal(LocalIndexType scvIdxLocal) const
    { return scvIndices_[scvIdxLocal]; }

    //! returns the local indices of the scvfs embedded in a local scv
    const LocalIndexContainer& localScvfIndicesInScv(LocalIndexType scvIdxLocal) const
    { return localScvfIndicesInScv_[scvIdxLocal]; }

    //! returns a global scvf idx for a given local scvf index
    IndexType scvfIdxGlobal(LocalIndexType scvfIdxLocal) const
    { return scvfIndices_[scvfIdxLocal]; }

    //! returns whether or not an scvf touches the boundary
    bool scvfIsOnBoundary(LocalIndexType scvfIdxLocal) const
    { return scvfIsOnBoundary_[scvfIdxLocal]; }

    //! returns the local indices of the neighboring scvs of an scvf
    const LocalIndexContainer& localNeighboringScvIndices(LocalIndexType scvfIdxLocal) const
    { return scvfNeighborScvIndices_[scvfIdxLocal]; }

private:
    //! returns the node-local scv index to a given global scv index
    unsigned int findLocalScvIdx_(IndexType globalScvIdx) const
    {
        auto it = std::find( scvIndices_.begin(), scvIndices_.end(), globalScvIdx );
        assert(it != scvIndices_.end() && "Global scv index not found in local container!");
        return std::distance(scvIndices_.begin(), it);
    }

    //! checks whether or not the inserted index set makes sense
    bool checkGlobalIndexSet_() const
    {
        //! do all scvs have dim embedded scvfs?
        for (const auto& scvfIndices : localScvfIndicesInScv_)
            if (scvfIndices.size() != dim)
                DUNE_THROW(Dune::InvalidStateException, "Number of scv faces found for an scv != dimension");

        //! are all scvfs unique?
        for (const auto& scvfIdx : scvfIndices_)
            for (const auto& scvfIdx2 : scvfIndices_)
                if (scvfIdx == scvfIdx2)
                    DUNE_THROW(Dune::InvalidStateException, "Inserted scvfs seem to not be unique");

        return true;
    }

    //! the indices of the scvs around a dual grid node
    GlobalIndexContainer scvIndices_;
    //! maps to each scv a list of scvf indices embedded in it
    std::vector<LocalIndexContainer> localScvfIndicesInScv_;

    //! the indices of the scvfs around a dual grid node
    GlobalIndexContainer scvfIndices_;
    //! maps to each scvf a boolean to indicate if it is on the boundary
    std::vector<bool> scvfIsOnBoundary_;
    //! maps to each scvf a list of neighbouring local scv indices
    //! ordering: 0 - inside scv idx; 1..n - outside scv indices
    std::vector<LocalIndexContainer> scvfNeighborScvIndices_;
    //! maps to each scvf a list of neighbouring global scv indices
    //! This container is destroyed when makeLocal() is called
    std::vector<GlobalIndexContainer> scvfNeighborScvIndicesGlobal_;
    //! stores how many boundary scvfs are embedded in this dual grid node
    std::size_t numBoundaryScvfs_;
};

/*!
 * \ingroup Mpfa
 * \brief Class for the index sets of the dual grid in mpfa schemes.
 */
template<class TypeTag>
class CCMpfaDualGridIndexSet : public std::vector<DualGridNodalIndexSet<TypeTag>>
{
    using ParentType = std::vector<DualGridNodalIndexSet<TypeTag>>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    static const int dim = GridView::dimension;

public:
    using NodalIndexSet = DualGridNodalIndexSet<TypeTag>;

    CCMpfaDualGridIndexSet(const GridView& gridView)
    {
        (*this).resize(gridView.size(dim));
    }

    //! after the global data has been inserted, this
    //! method should be called to set up the local index sets
    void makeLocalIndexSets()
    {
        for (auto& node : (*this))
            node.makeLocal();
    }

    //! access with an scvf
    const NodalIndexSet& operator[] (const SubControlVolumeFace& scvf) const
    { return (*this)[scvf.vertexIndex()]; }
    NodalIndexSet& operator[] (const SubControlVolumeFace& scvf)
    { return (*this)[scvf.vertexIndex()]; }

    //! access with an index
    using ParentType::operator[];
};
} // end namespace


#endif
