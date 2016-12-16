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
 * \brief Base class for the global interaction volumes of the mpfa-o method.
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_GLOBALINTERACTIONVOLUMESEEDS_BASE_HH
#define DUMUX_DISCRETIZATION_MPFA_GLOBALINTERACTIONVOLUMESEEDS_BASE_HH

#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include "globalfvgeometry.hh"

namespace Dumux
{
/*!
 * \ingroup MpfaO
 * \brief Base class for the creation and storage of the interaction volume seeds for mpfa methods.
 */
template<class TypeTag>
class CCMpfaGlobalInteractionVolumeSeedsBase
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, GlobalInteractionVolumeSeeds);
    // the actual implementation needs to be friend
    friend Implementation;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Helper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using InteractionVolumeSeed = typename InteractionVolume::Seed;
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using BoundaryInteractionVolumeSeed = typename BoundaryInteractionVolume::Seed;

    using IndexType = typename GridView::IndexSet::IndexType;

public:
    CCMpfaGlobalInteractionVolumeSeedsBase(const GridView gridView) : gridView_(gridView) {}

    // initializes the interaction volumes or the seeds
    void update(const Problem& problem, std::vector<bool>& boundaryVertices)
    {
        problemPtr_ = &problem;
        seeds_.clear();
        boundarySeeds_.clear();

        auto numScvf = problem_().model().globalFvGeometry().numScvf();
        scvfIndexMap_.resize(numScvf);
        std::vector<bool> isFaceHandled(numScvf, false);

        // initialize the seeds according to the mpfa method
        asImp_().initializeSeeds_(boundaryVertices, isFaceHandled);
    }

    const InteractionVolumeSeed& seed(const SubControlVolumeFace& scvf) const
    { return seeds_[scvfIndexMap_[scvf.index()]]; }

    const BoundaryInteractionVolumeSeed& boundarySeed(const SubControlVolumeFace& scvf) const
    { return boundarySeeds_[scvfIndexMap_[scvf.index()]]; }

private:

    void initializeSeeds_(std::vector<bool>& boundaryVertices,
                          std::vector<bool>& isFaceHandled)
    {
        const auto numScvf = problem_().model().globalFvGeometry().numScvf();
        const auto numBoundaryScvf = problem_().model().globalFvGeometry().numBoundaryScvf();

        // reserve memory
        seeds_.reserve(numScvf - numBoundaryScvf);
        boundarySeeds_.reserve(numBoundaryScvf);

        IndexType boundarySeedIndex = 0;
        IndexType seedIndex = 0;
        for (const auto& element : elements(gridView_))
        {
            auto fvGeometry = localView(problem_().model().globalFvGeometry());
            fvGeometry.bindElement(element);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                // skip the rest if we already handled this face
                if (isFaceHandled[scvf.index()])
                    continue;

                if (boundaryVertices[scvf.vertexIndex()])
                {
                    // make the boundary interaction volume seed
                    boundarySeeds_.emplace_back(Helper::makeBoundaryInteractionVolumeSeed(problem_(), element, fvGeometry, scvf));

                    // update the index map entries for the global scv faces in the interaction volume
                    for (auto scvfIdxGlobal : boundarySeeds_.back().globalScvfIndices())
                    {
                        scvfIndexMap_[scvfIdxGlobal] = boundarySeedIndex;
                        isFaceHandled[scvfIdxGlobal] = true;
                    }

                    // increment counter
                    boundarySeedIndex++;
                }
                else
                {
                    // make the inner interaction volume seed
                    seeds_.emplace_back(Helper::makeInnerInteractionVolumeSeed(problem_(), element, fvGeometry, scvf));

                    // update the index map entries for the global scv faces in the interaction volume
                    for (auto scvfIdxGlobal : seeds_.back().globalScvfIndices())
                    {
                        scvfIndexMap_[scvfIdxGlobal] = seedIndex;
                        isFaceHandled[scvfIdxGlobal] = true;
                    }

                    // increment counter
                    seedIndex++;
                }
            }
        }

        // shrink vectors to actual size
        seeds_.shrink_to_fit();
        boundarySeeds_.shrink_to_fit();
    }

    const Problem& problem_() const
    { return *problemPtr_; }

    const Implementation& asImp_() const
    { return *static_cast<Implementation*>(this); }

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Problem* problemPtr_;
    GridView gridView_;
    std::vector<IndexType> scvfIndexMap_;
    std::vector<InteractionVolumeSeed> seeds_;
    std::vector<BoundaryInteractionVolumeSeed> boundarySeeds_;
};
} // end namespace


#endif
