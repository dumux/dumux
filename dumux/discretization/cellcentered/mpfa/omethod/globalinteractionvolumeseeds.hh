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
#ifndef DUMUX_DISCRETIZATION_MPFA_O_GLOBALINTERACTIONVOLUMESEEDS_HH
#define DUMUX_DISCRETIZATION_MPFA_O_GLOBALINTERACTIONVOLUMESEEDS_HH

#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>

namespace Dumux
{
/*!
 * \ingroup MpfaO
 * \brief Base class for the creation and storage of the interaction volume seeds for mpfa methods.
 */
template<class TypeTag>
class CCMpfaGlobalInteractionVolumeSeedsImplementation<TypeTag, MpfaMethods::oMethod>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Helper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using InteractionVolumeSeed = typename InteractionVolume::Seed;
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using BoundaryInteractionVolumeSeed = typename BoundaryInteractionVolume::Seed;

    using IndexType = typename GridView::IndexSet::IndexType;
    using LocalIndexType = typename InteractionVolume::LocalIndexType;

    static const int dim = GridView::dimension;
    static const int numEq = GET_PROP_VALUE(TypeTag, NumEq);

public:
    CCMpfaGlobalInteractionVolumeSeedsImplementation(const GridView gridView) : gridView_(gridView) {}

    // initializes the interaction volumes or the seeds
    void update(const Problem& problem)
    {
        problemPtr_ = &problem;
        seeds_.clear();
        boundarySeeds_.clear();
        scvfIndexMap_.clear();
        initializeSeeds_();
    }

    const InteractionVolumeSeed& seed(const SubControlVolumeFace& scvf) const
    {
        return seeds_[scvfIndexMap_[scv.index()]];
    }

    const BoundaryInteractionVolume& boundarySeed(const SubControlVolumeFace& scvf, const LocalIndexType eqIdx) const
    {
        assert(boundarySeeds_[scvfIndexMap_[scv.index()]].size() > eqIdx);
        return boundarySeeds_[scvfIndexMap_[scv.index()]][eqIdx];
    }

private:
    void initializeSeeds_()
    {
        seeds_.reserve(gridView_.size(dim));

        auto numScvf = problem_().model().fvGeometries().numScvf();
        scvfIndexMap_.resize(numScvf, -1);
        scvfTouchesBoundary_.resize(numScvf, false);

        // first we have to construct all the interaction volumes that touch an inner or outer boundary
        IndexType boundarySeedIndex = 0;
        for (const auto& element : element(gridView_))
        {
            auto fvGeometry = localView(problem_().model().globalFvGeometries());

            // at this point we can only use the bind method, as the stencils are not yet set up
            fvGeometry.bind(element);
            for (const auto& scvf : fvGeometry.scvfs())
            {
                // skip the rest if we already handled this face
                if (scvfIndexMap_[scvf.index()] != -1)
                    continue;

                // also, skip the rest if this face doesn't touch a boundary
                bool touchesBoundary = false;
                for (LocalIndexType eqIdx = 0; eqIdx < numEq; ++eqIdx)
                {
                    if (Helper::getMpfaFaceType(problem_(), element, scvf, eqIdx) != MpfaFaceTypes::interior)
                    {
                        touchesBoundary = true;
                        break;
                    }
                }

                if (!touchesBoundary)
                    continue;

                // container to store the interaction volume seeds
                std::vector<BoundaryInteractionVolume> seedVector;
                seedVector.reserve(numEq);
                for (LocalIndexType eqIdx = 0; eqIdx < numEq; ++eqIdx)
                    seedVector.emplace_back(Helper::makeBoundaryInteractionVolumeSeed(problem_(), element, fvGeometry, scvf, eqIdx));

                // update the index map entries for the global scv faces in the interaction volume
                for (const auto& localScvf : seedVector[0].scvfSeeds())
                    for (const auto scvfIdxGlobal : localScvf.globalScvfIndices())
                    {
                        assert(scvfIndexMap_[scvfIdxGlobal] == -1);
                        scvfIndexMap_[scvfIdxGlobal] = boundarySeedIndex;
                        scvfTouchesBoundary_[scvfIdxGlobal] = true;
                    }

                // store interaction volume and increment counter
                boundarySeeds_.emplace_back(std::move(seedVector));
                boundarySeedIndex++;
            }
        }

        // Now we construct the inner interaction volume seeds
        IndexType seedIndex = 0;
        for (const auto& element : element(gridView_))
        {
            auto fvGeometry = localView(problem_().model().globalFvGeometries());

            // at this point we can only use the bind method, as the stencils are not yet set up
            fvGeometry.bind(element);
            for (const auto& scvf : fvGeometry.scvfs())
            {
                if (scvfIndexMap_[scvf.index()] != -1)
                    continue;

                // the inner interaction volume seed
                auto seed = Helper::makeInnerInteractionVolumeSeed(problem_(), element, fvGeometry, scvf);

                // update the index map entries for the global scv faces in the interaction volume
                for (const auto& localScvf : seedVector[0].scvfSeeds())
                    for (const auto scvfIdxGlobal : localScvf.globalScvfIndices())
                    {
                        assert(scvfIndexMap_[scvfIdxGlobal] == -1);
                        scvfIndexMap_[scvfIdxGlobal] = seedIndex;
                    }

                // store interaction volume and increment counter
                seeds_.emplace_back(std::move(seed));
                seedIndex++;
            }
        }
    }

    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;
    GridView gridView_;
    std::vector<IndexType> scvfIndexMap_;
    std::vector<bool> scvfTouchesBoundary_;
    std::vector<InteractionVolumeSeed> seeds_;
    std::vector<std::vector<BoundaryInteractionVolumeSeed>> boundarySeeds_;
};

} // end namespace

#endif
