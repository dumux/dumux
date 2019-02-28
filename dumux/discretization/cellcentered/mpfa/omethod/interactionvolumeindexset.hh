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
 * \ingroup CCMpfaDiscretization
 * \brief Class for the index set within an interaction volume of the mpfa-o scheme.
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_O_INTERACTIONVOLUME_INDEXSET_HH
#define DUMUX_DISCRETIZATION_MPFA_O_INTERACTIONVOLUME_INDEXSET_HH

#include <dune/common/reservedvector.hh>

#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief The interaction volume index set class for the mpfa-o scheme.
 *
 * \tparam DualGridNodalIndexSet The type used for the nodal index set in the dual grid.
 */
template< class DualGridNodalIndexSet >
class CCMpfaOInteractionVolumeIndexSet
{
public:
    //! Export the type used for the nodal grid index sets
    using NodalIndexSet = DualGridNodalIndexSet;

    //! Export the types used for local/grid indices
    using LocalIndexType = typename DualGridNodalIndexSet::LocalIndexType;
    using GridIndexType = typename DualGridNodalIndexSet::GridIndexType;

    //! Export the stencil types used
    using NodalGridStencilType = typename DualGridNodalIndexSet::NodalGridStencilType;
    using NodalLocalStencilType = typename DualGridNodalIndexSet::NodalLocalStencilType;
    using NodalGridScvfStencilType = typename DualGridNodalIndexSet::NodalGridScvfStencilType;

    //! Export the type used for the neighbor scv index sets of the scvfs
    using ScvfNeighborLocalIndexSet = typename DualGridNodalIndexSet::ScvfNeighborLocalIndexSet;

    //! The constructor
    template< class FlipScvfIndexSet >
    CCMpfaOInteractionVolumeIndexSet(const NodalIndexSet& nodalIndexSet, const FlipScvfIndexSet& flipScvfIndexSet)
    : nodalIndexSet_(nodalIndexSet)
    {
        const auto numNodalScvfs = nodalIndexSet.numScvfs();

        // kee track of which nodal scvfs have been handled already
        nodeToIvScvf_.resize(numNodalScvfs);
        std::vector<bool> isHandled(numNodalScvfs, false);

        // go over faces in nodal index set, check if iv-local face has been
        // inserted already for this scvf and if not, insert index mapping
        numFaces_ = 0;
        for (LocalIndexType i = 0; i < numNodalScvfs; ++i)
        {
            // check if the nodal scvf still has to be handled
            if (isHandled[i])
                continue;

            // for scvfs touching the boundary there are no "outside" scvfs
            if (nodalIndexSet.scvfIsOnBoundary(i))
            {
                scvfNeighborScvLocalIndices_.push_back({nodalIndexSet.insideScvLocalIndex(i)});
                nodeToIvScvf_[i] = ivToNodeScvf_.size();
                isHandled[i] = true;
                ivToNodeScvf_.push_back(i);
                numFaces_++;
                continue;
            }

            // insert a new iv-local face
            const auto curIvLocalScvfIdx = ivToNodeScvf_.size();
            nodeToIvScvf_[i] = curIvLocalScvfIdx;
            isHandled[i] = true;

            // construct local index sets
            const auto& flipScvfIndices = flipScvfIndexSet[nodalIndexSet.gridScvfIndex(i)];
            const auto numFlipIndices = flipScvfIndices.size();

            ScvfNeighborLocalIndexSet neighborsLocal;
            neighborsLocal.resize(numFlipIndices + 1);
            neighborsLocal[0] = nodalIndexSet.insideScvLocalIndex(i);

            // mappings for all flip scvf
            for (unsigned int j = 0; j < numFlipIndices; ++j)
            {
                const auto outsideScvfIdx = flipScvfIndices[j];
                for (unsigned int nodeLocalIdx = 0; nodeLocalIdx < nodalIndexSet.numScvfs(); ++nodeLocalIdx)
                {
                    if (nodalIndexSet.gridScvfIndex(nodeLocalIdx) == outsideScvfIdx)
                    {
                        neighborsLocal[j+1] = nodalIndexSet.insideScvLocalIndex(nodeLocalIdx);
                        nodeToIvScvf_[nodeLocalIdx] = curIvLocalScvfIdx;
                        isHandled[nodeLocalIdx] = true;
                        break; // go to next outside scvf
                    }
                }
            }

            scvfNeighborScvLocalIndices_.push_back(neighborsLocal);
            ivToNodeScvf_.push_back(i);
            numFaces_++;
        }
    }

    //! returns the corresponding nodal index set
    const NodalIndexSet& nodalIndexSet() const
    { return nodalIndexSet_; }

    //! returns the global scv indices connected to this dual grid node
    const NodalGridStencilType& gridScvIndices() const
    { return nodalIndexSet_.gridScvIndices(); }

    //! returns the global scvf indices embedded in this interaction volume
    const NodalGridScvfStencilType& gridScvfIndices() const
    { return nodalIndexSet_.gridScvfIndices(); }

    //! returns the number of faces in the interaction volume
    std::size_t numFaces() const
    { return numFaces_; }

    //! returns the number of scvs in the interaction volume
    std::size_t numScvs() const
    { return nodalIndexSet_.numScvs(); }

    //! returns a grid scv idx for a given iv-local scv index
    GridIndexType gridScvIndex(LocalIndexType ivLocalScvIdx) const
    {
        assert(ivLocalScvIdx < numScvs());
        return gridScvIndices()[ivLocalScvIdx];
    }

    //! returns a grid scvf idx for a given iv-local scvf index
    GridIndexType gridScvfIndex(LocalIndexType ivLocalScvfIdx) const
    {
        assert(ivLocalScvfIdx < numFaces());
        return nodalIndexSet_.gridScvfIndex( ivToNodeScvf_[ivLocalScvfIdx] );
    }

    //! returns the iv-local scvf idx of the i-th scvf embedded in a local scv
    LocalIndexType localScvfIndex(LocalIndexType scvIdxLocal, unsigned int i) const
    {
        assert(nodalIndexSet_.localScvfIndex(scvIdxLocal, i) < nodeToIvScvf_.size());
        return nodeToIvScvf_[ nodalIndexSet_.localScvfIndex(scvIdxLocal, i) ];
    }

    //! returns the local indices of the neighboring scvs of an scvf
    const ScvfNeighborLocalIndexSet& neighboringLocalScvIndices(LocalIndexType ivLocalScvfIdx) const
    {
        assert(ivLocalScvfIdx < numFaces());
        return scvfNeighborScvLocalIndices_[ivLocalScvfIdx];
    }

private:
    using NI = NodalIndexSet;

    std::size_t numFaces_;
    const NI& nodalIndexSet_;
    // Index maps from and to nodal index set. For the map to the
    // nodal set we use the same storage type as we know the nodal
    // has more faces, thus sufficient guaranteed here!
    typename NI::Traits::template NodalScvfDataStorage< LocalIndexType > ivToNodeScvf_;
    typename NI::Traits::template NodalScvfDataStorage< LocalIndexType > nodeToIvScvf_;
    // maps to each scvf a list of neighbouring scv indices
    // ordering: 0 - inside scv idx; 1..n - outside scv indices
    typename NI::Traits::template NodalScvfDataStorage< ScvfNeighborLocalIndexSet > scvfNeighborScvLocalIndices_;
};

} // end namespace Dumux

#endif
