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
 * \brief Class for the index set within an interaction volume of the mpfa-o scheme.
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_O_NTERACTIONVOLUME_INDEXSET_HH
#define DUMUX_DISCRETIZATION_MPFA_O_INTERACTIONVOLUME_INDEXSET_HH

#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>

namespace Dumux
{
/*!
 * \ingroup Mpfa
 * \brief The interaction volume index set class for the mpfa-o scheme.
 */
template<class DualGridNodalIndexSet, class GlobalIndexContainer, class LocalIndexContainer>
class CCMpfaOInteractionVolumeIndexSet
{
    using LocalIndexType = typename LocalIndexContainer::value_type;
    using GlobalIndexType = typename GlobalIndexContainer::value_type;

public:
    CCMpfaOInteractionVolumeIndexSet(const DualGridNodalIndexSet& nodalIndexSet) : nodalIndexSet_(nodalIndexSet)
    {
        //! determine the number of iv-local faces
        const auto numNodalScvfs = nodalIndexSet.numScvfs();
        const auto numBoundaryScvfs = nodalIndexSet.numBoundaryScvfs();
        numFaces_ = numBoundaryScvfs + (numNodalScvfs-numBoundaryScvfs)/2;

        // make sure we found a reasonable number of faces
        assert((numNodalScvfs-numBoundaryScvfs)%2 == 0);

        //! index transformation from interaction-volume-local to node-local
        ivToNodeScvf_.reserve(numFaces_);
        nodeToIvScvf_.resize(numNodalScvfs);

        //! keeps track of which nodal scvfs have been handled already
        std::vector<bool> isHandled(numNodalScvfs, false);

        //! go over faces in nodal index set, check if iv-local face
        //! has been inserted already for this scvf and if not, insert
        //! the index transformation
        unsigned int findCounter = 0;
        for (LocalIndexType i = 0; i < numNodalScvfs; ++i)
        {
            //! check if the nodal scvf still has to be handled
            if (isHandled[i])
                continue;

            //! for scvfs touching the boundary there are no "outside" scvfs
            if (nodalIndexSet.scvfIsOnBoundary(i))
            {
                nodeToIvScvf_[i] = ivToNodeScvf_.size();
                isHandled[i] = true;
                ivToNodeScvf_.push_back(i);
                findCounter++; continue;
            }

            //! We insert a new iv-local face and find all "outside" scvfs that map
            //! to this face as well by comparing the set of neighboring scv indices.
            const auto scvIndices = [&nodalIndexSet, i] ()
                                    {
                                        auto tmp = nodalIndexSet.localNeighboringScvIndices(i);
                                        std::sort(tmp.begin(), tmp.end());
                                        return tmp;
                                    } ();
            const auto numNeighborsI = scvIndices.size();

            std::vector<LocalIndexType> outsideScvfs;
            for (LocalIndexType j = i+1; j < numNodalScvfs; ++j)
            {
                //! a face that has been handled already cannot be an "outside" face here
                if (!isHandled[j])
                {
                    const auto scvIndices2 = [&nodalIndexSet, j] ()
                                             {
                                                 auto tmp = nodalIndexSet.localNeighboringScvIndices(j);
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
            findCounter++;
        }

        //! ensure we found as many faces as anticipated
        assert(findCounter == numFaces_ && "Couldn't find as many faces as anticipated");
    }

    //! returns the corresponding nodal index set
    const DualGridNodalIndexSet& nodalIndexSet() const
    { return nodalIndexSet_; }

    //! returns a global scvf idx for a given iv_local scvf index
    GlobalIndexType scvfIdxGlobal(LocalIndexType ivLocalScvfIdx) const
    { return nodalIndexSet_.scvfIdxGlobal(ivToNodeScvf_[ivLocalScvfIdx]); }

    //! returns the iv-local scvf idx of the i-th (coordDir) scvfs embedded in a local scv
    LocalIndexType localScvfIndexInScv(LocalIndexType scvIdxLocal, unsigned int coordDir) const
    { return nodeToIvScvf_[nodalIndexSet_.localScvfIndicesInScv(scvIdxLocal)[coordDir]]; }

    //! returns whether or not an scvf touches the boundary
    bool scvfTouchesBoundary(LocalIndexType ivLocalScvfIdx) const
    { return nodalIndexSet_.scvfTouchesBoundary(ivToNodeScvf_[ivLocalScvfIdx]); }

    //! returns the local indices of the neighboring scvs of an scvf
    const LocalIndexContainer& localNeighboringScvIndices(LocalIndexType ivLocalScvfIdx) const
    { return nodalIndexSet_.localNeighboringScvIndices(ivToNodeScvf_[ivLocalScvfIdx]); }

    //! returns the number of faces in the interaction volume
    std::size_t numFaces() const
    { return numFaces_; }

    //! returns the number of scvs in the interaction volume
    std::size_t numScvs() const
    { return nodalIndexSet_.numScvs(); }

private:
    const DualGridNodalIndexSet& nodalIndexSet_;

    std::size_t numFaces_;
    LocalIndexContainer ivToNodeScvf_;
    LocalIndexContainer nodeToIvScvf_;
};
} // end namespace


#endif
