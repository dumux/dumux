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
 * \ingroup CCMpfaDiscretization
 * \brief Class for the index set within an interaction volume of the mpfa-o scheme.
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_O_INTERACTIONVOLUME_INDEXSET_HH
#define DUMUX_DISCRETIZATION_MPFA_O_INTERACTIONVOLUME_INDEXSET_HH

#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>

namespace Dumux
{
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
    using LocalIndexType = typename DualGridNodalIndexSet::LocalIndexType;
    using GridIndexType = typename DualGridNodalIndexSet::GridIndexType;

    using LocalIndexContainer = typename DualGridNodalIndexSet::LocalIndexContainer;
    using GridIndexContainer = typename DualGridNodalIndexSet::GridIndexContainer;

    /*!
     * \brief The constructor
     * \note The actual type used for the nodal index sets might be different, as maybe
     *       a different type for the local indexes is used. We therefore template this
     *       constructor. However, a static assertion enforces you to use the same LocalIndexType
     *       in the traits for both the secondary and the primary interaction volume traits.
     *
     * \tparam NodalIndexSet Possibly differing type for the DualGridNodalIndexSet
     */
    template< class NodalIndexSet >
    CCMpfaOInteractionVolumeIndexSet(const NodalIndexSet& nodalIndexSet)
    : nodalIndexSet_( static_cast<const DualGridNodalIndexSet&>(nodalIndexSet) )
    {
        // make sure the index types used are the same in order to avoid any losses due to type conversion
        static_assert(std::is_same<GridIndexType, typename NodalIndexSet::GridIndexType>::value,
                      "Provided nodal index set does not use the same type for grid indices as the given template argument");
        static_assert(std::is_same<LocalIndexType, typename NodalIndexSet::LocalIndexType>::value,
                      "Provided nodal index set does not use the same type for local indices as the given template argument");

        // determine the number of iv-local faces for memory reservation
        // note that this might be a vast overestimation on surface grids!
        const auto numNodalScvfs = nodalIndexSet.numScvfs();
        const auto numBoundaryScvfs = nodalIndexSet.numBoundaryScvfs();
        const std::size_t numFaceEstimate = numBoundaryScvfs + (numNodalScvfs-numBoundaryScvfs)/2;

        // make sure we found a reasonable number of faces
        assert((numNodalScvfs-numBoundaryScvfs)%2 == 0);

        // index transformation from interaction-volume-local to node-local
        ivToNodeScvf_.reserve(numFaceEstimate);
        nodeToIvScvf_.resize(numNodalScvfs);

        // the local neighboring scv indices of the faces
        scvfNeighborScvLocalIndices_.reserve(numFaceEstimate);

        // keeps track of which nodal scvfs have been handled already
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
                nodeToIvScvf_[i] = ivToNodeScvf_.size();
                isHandled[i] = true;
                ivToNodeScvf_.push_back(i);
                numFaces_++;
                continue;
            }

            // We insert a new iv-local face and find all "outside" scvfs that map
            // to this face as well by comparing the set of neighboring scv indices.
            const auto scvIndices = [&nodalIndexSet, i] ()
                                    {
                                        auto tmp = nodalIndexSet.neighboringScvIndices(i);
                                        std::sort(tmp.begin(), tmp.end());
                                        return tmp;
                                    } ();
            const auto numNeighborsI = scvIndices.size();

            std::vector<LocalIndexType> outsideScvfs;
            for (LocalIndexType j = i+1; j < numNodalScvfs; ++j)
            {
                // a face that has been handled already cannot be an "outside" face here
                if (!isHandled[j])
                {
                    const auto scvIndices2 = [&nodalIndexSet, j] ()
                                             {
                                                 auto tmp = nodalIndexSet.neighboringScvIndices(j);
                                                 std::sort(tmp.begin(), tmp.end());
                                                 return tmp;
                                             } ();

                    // if the sizes aren't equal, this can't be an "outside" face
                    if (scvIndices2.size() != numNeighborsI)
                        continue;

                    // if the neighboring scv indices are the same, this is an "outside" face
                    if (std::equal(scvIndices.begin(), scvIndices.end(), scvIndices2.begin()))
                        outsideScvfs.push_back(j);
                }
            }

            // insert mappings
            const auto curIvLocalScvfIdx = ivToNodeScvf_.size();
            nodeToIvScvf_[i] = curIvLocalScvfIdx;
            isHandled[i] = true;
            for (const auto nodeLocalScvfIdx : outsideScvfs)
            {
                nodeToIvScvf_[nodeLocalScvfIdx] = curIvLocalScvfIdx;
                isHandled[nodeLocalScvfIdx] = true;
            }
            ivToNodeScvf_.push_back(i);
            numFaces_++;
        }

        // compute local neighboring scv indices for the iv-local scvfs
        scvfNeighborScvLocalIndices_.resize(numFaces_);
        for (unsigned int i = 0; i < numFaces_; ++i)
        {
            const auto& neighborsGlobal = nodalIndexSet_.neighboringScvIndices(ivToNodeScvf_[i]);
            const auto numNeighbors = nodalIndexSet_.scvfIsOnBoundary(ivToNodeScvf_[i]) ? 1 : neighborsGlobal.size();

            scvfNeighborScvLocalIndices_[i].resize(numNeighbors);
            for (unsigned int j = 0; j < numNeighbors; ++j)
                scvfNeighborScvLocalIndices_[i][j] = findLocalScvIdx_(neighborsGlobal[j]);
        }
    }

    //! returns the corresponding nodal index set
    const DualGridNodalIndexSet& nodalIndexSet() const { return nodalIndexSet_; }

    //! returns a global scvf idx for a given iv_local scvf index
    GridIndexType scvfIdxGlobal(LocalIndexType ivLocalScvfIdx) const
    { return nodalIndexSet_.scvfIdxGlobal( ivToNodeScvf_[ivLocalScvfIdx] ); }

    //! returns the iv-local scvf idx of the i-th scvfs embedded in a local scv
    LocalIndexType scvfIdxLocal(LocalIndexType scvIdxLocal, unsigned int i) const
    { return nodeToIvScvf_[ nodalIndexSet_.scvfIdxLocal(scvIdxLocal, i) ]; }

    //! returns the local indices of the neighboring scvs of an scvf
    const LocalIndexContainer& neighboringLocalScvIndices(LocalIndexType ivLocalScvfIdx) const
    { return scvfNeighborScvLocalIndices_[ivLocalScvfIdx]; }

    //! returns the number of faces in the interaction volume
    std::size_t numFaces() const { return numFaces_; }

    //! returns the number of scvs in the interaction volume
    std::size_t numScvs() const { return nodalIndexSet_.numScvs(); }

    //! returns the global scv indices connected to this dual grid node
    const GridIndexContainer& globalScvIndices() const { return nodalIndexSet_.globalScvIndices(); }

    //! returns the global scvf indices connected to this dual grid node
    const GridIndexContainer& globalScvfIndices() const { return nodalIndexSet_.globalScvfIndices(); }

private:
    //! returns the local scv index to a given global scv index
    unsigned int findLocalScvIdx_(GridIndexType globalScvIdx) const
    {
        auto it = std::find( nodalIndexSet_.globalScvIndices().begin(), nodalIndexSet_.globalScvIndices().end(), globalScvIdx );
        assert(it != nodalIndexSet_.globalScvIndices().end() && "Global scv index not found in local container!");
        return std::distance(nodalIndexSet_.globalScvIndices().begin(), it);
    }

    const DualGridNodalIndexSet& nodalIndexSet_;

    std::size_t numFaces_;
    std::vector<LocalIndexType> ivToNodeScvf_;
    std::vector<LocalIndexType> nodeToIvScvf_;
    // maps to each scvf a list of neighbouring scv indices
    // ordering: 0 - inside scv idx; 1..n - outside scv indices
    std::vector< LocalIndexContainer > scvfNeighborScvLocalIndices_;
};

} // end namespace Dumux

#endif
