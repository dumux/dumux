// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCMpfaDiscretization
 * \brief Index set for the dual grid in mpfa schemes.
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_DUALGRID_INDEX_SET_HH
#define DUMUX_DISCRETIZATION_MPFA_DUALGRID_INDEX_SET_HH

#include <cassert>
#include <vector>
#include <algorithm>
#include <cstdint>

#include <dumux/common/indextraits.hh>
#include <dumux/common/reservedvector.hh>

namespace Dumux {
namespace CCMpfa {

template<typename GridView>
inline constexpr unsigned int minNumScvsAtVertex = int(GridView::dimension) == 2 ? 4 : 8;

template<typename GridView>
inline constexpr unsigned int minNumScvfsAtVertex = int(GridView::dimension) == 2 ? 4 : 12;

template<typename GridView>
struct DataStorage
{
    template<typename T> using NodalScvDataStorage = Dumux::ReservedVector<T, minNumScvsAtVertex<GridView>>;
    template<typename T> using NodalScvfDataStorage = Dumux::ReservedVector<T, minNumScvfsAtVertex<GridView>>;
    template<typename T> using ScvfNeighborDataStorage = std::conditional_t<
        int(GridView::dimension) < int(GridView::dimensionworld),
        std::vector<T>,
        Dumux::ReservedVector<T, 2>
    >;
};

} // namespace Mpfa

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Nodal index set for mpfa schemes, constructed around grid vertices.
 * \tparam GV The grid view
 */
template<class GV>
class CCMpfaDualGridNodalIndexSet
{
    using Storage = CCMpfa::DataStorage<GV>;

    using LI = std::uint_least16_t;
    using GI = typename GV::IndexSet::IndexType;
    using DimIndexVector = Dumux::ReservedVector<LI, GV::dimension>;
    using ScvfIndicesInScvStorage = typename Storage::NodalScvDataStorage<DimIndexVector>;

public:
    using LocalIndexType = LI;
    using GridIndexType = GI;

    // for compatibility
    struct Traits {
        using GridView = GV;
    };

    //! Constructor
    CCMpfaDualGridNodalIndexSet() : numBoundaryScvfs_(0) {}

    //! Inserts data for a given scvf
    template<typename SubControlVolumeFace>
    void insert(const SubControlVolumeFace& scvf)
    { insert(scvf.index(), scvf.insideScvIdx(), scvf.boundary()); }

    //! Inserts scvf data
    void insert(const GridIndexType scvfIdx,
                const GridIndexType insideScvIdx,
                const bool boundary)
    {
        // this should always be called only once per scvf
        assert(std::count(scvfIndices_.begin(), scvfIndices_.end(), scvfIdx ) == 0 && "scvf already inserted!");

        // the local index of the scvf data about to be inserted
        const LocalIndexType curScvfLocalIdx = scvfIndices_.size();

        // if scvf is on boundary, increase counter
        if (boundary) numBoundaryScvfs_++;

        // insert data on the new scv
        scvfIndices_.push_back(scvfIdx);
        scvfIsOnBoundary_.push_back(boundary);

        // if entry for the inside scv exists append data, create otherwise
        auto it = std::find( scvIndices_.begin(), scvIndices_.end(), insideScvIdx );
        if (it != scvIndices_.end())
        {
            const auto scvIdxLocal = std::distance(scvIndices_.begin(), it);
            scvfInsideScvIndices_.push_back(scvIdxLocal);
            localScvfIndicesInScv_[scvIdxLocal].push_back(curScvfLocalIdx);
        }
        else
        {
            scvfInsideScvIndices_.push_back(scvIndices_.size());
            localScvfIndicesInScv_.push_back({curScvfLocalIdx});
            scvIndices_.push_back(insideScvIdx);
        }
    }

    //! returns the number of scvs around the node
    std::size_t numScvs() const
    { return scvIndices_.size(); }

    //! returns the number of scvfs around the node
    std::size_t numScvfs() const
    { return scvfIndices_.size(); }

    //! returns the number of boundary scvfs around the node
    std::size_t numBoundaryScvfs() const
    { return numBoundaryScvfs_; }

    //! returns the grid scv indices connected to this dual grid node
    const auto& gridScvIndices() const
    { return scvIndices_; }

    //! returns the grid scvf indices connected to this dual grid node
    const auto& gridScvfIndices() const
    { return scvfIndices_; }

    //! returns whether or not the i-th scvf is on a domain boundary
    bool scvfIsOnBoundary(unsigned int i) const
    {
        assert(i < numScvfs());
        return scvfIsOnBoundary_[i];
    }

    //! returns the grid scv idx of the i-th scv
    GridIndexType gridScvIndex(unsigned int i) const
    {
        assert(i < numScvs());
        return scvIndices_[i];
    }

    //! returns the index of the i-th scvf
    GridIndexType gridScvfIndex(unsigned int i) const
    {
        assert(i < numScvfs());
        return scvfIndices_[i];
    }

    //! returns the grid index of the j-th scvf embedded in the i-th scv
    GridIndexType gridScvfIndex(unsigned int i, unsigned int j) const
    {
        assert(i < numScvs());
        assert(j < localScvfIndicesInScv_[i].size());
        return scvfIndices_[ localScvfIndicesInScv_[i][j] ];
    }

    //! returns the node-local index of the j-th scvf embedded in the i-th scv
    LocalIndexType localScvfIndex(unsigned int i, unsigned int j) const
    {
        assert(i < numScvs());
        assert(j < localScvfIndicesInScv_[i].size());
        return localScvfIndicesInScv_[i][j];
    }

    //! returns the node-local index of the inside scv of the i-th scvf
    LocalIndexType insideScvLocalIndex(unsigned int i) const
    {
        assert(i < numScvfs());
        return scvfInsideScvIndices_[i];
    }

private:
    typename Storage::NodalScvDataStorage<GI> scvIndices_;                        //!< The indices of the scvs around a dual grid node
    typename Storage::NodalScvDataStorage<DimIndexVector> localScvfIndicesInScv_; //!< Maps to each scv a list of scvf indices embedded in it

    std::size_t numBoundaryScvfs_;                                    //!< stores how many boundary scvfs are embedded in this dual grid node
    typename Storage::NodalScvfDataStorage<GI> scvfIndices_;          //!< the indices of the scvfs around a dual grid node
    typename Storage::NodalScvfDataStorage<bool> scvfIsOnBoundary_;   //!< Maps to each scvf a boolean to indicate if it is on the boundary
    typename Storage::NodalScvfDataStorage<LI> scvfInsideScvIndices_; //!< The inside local scv index for each scvf
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Class for the index sets of the dual grid in mpfa schemes.
 */
template<class GridView>
class CCMpfaDualGridIndexSet
{
public:
    using NodalIndexSet = CCMpfaDualGridNodalIndexSet<GridView>;
    using GridIndexType = typename NodalIndexSet::GridIndexType;

    CCMpfaDualGridIndexSet() = delete;
    CCMpfaDualGridIndexSet(const GridView& gridView)
    : nodalIndexSets_(gridView.size(GridView::dimension))
    {}

    template<class SubControlVolumeFace>
    const NodalIndexSet& operator[](const SubControlVolumeFace& scvf) const
    { return nodalIndexSets_[scvf.vertexIndex()]; }

    template<class SubControlVolumeFace>
    NodalIndexSet& operator[](const SubControlVolumeFace& scvf)
    { return nodalIndexSets_[scvf.vertexIndex()]; }

    const NodalIndexSet& operator[](GridIndexType i) const
    { return nodalIndexSets_[i]; }

    NodalIndexSet& operator[](GridIndexType i)
    { return nodalIndexSets_[i]; }

private:
    std::vector<NodalIndexSet> nodalIndexSets_;
};

} // end namespace Dumux

#endif
