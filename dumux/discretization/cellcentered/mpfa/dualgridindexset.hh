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
        assert( std::find(scvfIndices_.begin(), scvfIndices_.end(), scvfIdx ) == scvfIndices_.end() && "scvf has already been inserted!");

        // the local index of the scvf data about to be inserted
        const LocalIndexType curScvfLocalIdx = scvfIndices_.size();

        // add global scvf data
        GlobalIndexContainer scvIndices;
        scvIndices.reserve(outsideScvIndices.size()+1);
        scvIndices.push_back(insideScvIdx);
        scvIndices.insert(scvIndices.end(), outsideScvIndices.begin(), outsideScvIndices.end());

        // if scvf is on boundary, increase counter
        if (boundary) numBoundaryScvfs_++;

        // insert data on the new scv
        scvfIndices_.push_back(scvfIdx);
        scvfIsOnBoundary_.push_back(boundary);
        scvfNeighborScvIndices_.emplace_back( std::move(scvIndices) );

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

    //! returns the global scv idx of the i-th scv
    IndexType scvIdxGlobal(unsigned int i) const
    { return scvIndices_[i]; }

    //! returns the global index of the j-th scvf embedded in the i-th scv
    IndexType scvfIdxGlobal(unsigned int i, unsigned int j) const
    {
        assert(j < dim); // only dim faces can be embedded in an scv
        return scvfIndices_[ localScvfIndicesInScv_[i][j] ];
    }

    //! returns the nodel-local index of the j-th scvf embedded in the i-th scv
    IndexType scvfIdxLocal(unsigned int i, unsigned int j) const
    {
        assert(j < dim); // only dim faces can be embedded in an scv
        return localScvfIndicesInScv_[i][j];
    }

    //! returns the index of the i-th scvf
    IndexType scvfIdxGlobal(unsigned int i) const
    { return scvfIndices_[i]; }

    //! returns whether or not the i-th scvf touches the boundary
    bool scvfIsOnBoundary(unsigned int i) const
    { return scvfIsOnBoundary_[i]; }

    //! returns the indices of the neighboring scvs of the i-th scvf
    const GlobalIndexContainer& neighboringScvIndices(unsigned int i) const
    { return scvfNeighborScvIndices_[i]; }

private:
    //! the indices of the scvs around a dual grid node
    GlobalIndexContainer scvIndices_;
    //! maps to each scv a list of scvf indices embedded in it
    std::vector<LocalIndexContainer> localScvfIndicesInScv_;

    //! the indices of the scvfs around a dual grid node
    GlobalIndexContainer scvfIndices_;
    //! maps to each scvf a boolean to indicate if it is on the boundary
    std::vector<bool> scvfIsOnBoundary_;
    //! maps to each scvf a list of neighbouring scv indices
    //! ordering: 0 - inside scv idx; 1..n - outside scv indices
    std::vector<GlobalIndexContainer> scvfNeighborScvIndices_;
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
