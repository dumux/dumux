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
#ifndef DUMUX_DISCRETIZATION_MPFA_O_HYBRIDFPS_GLOBALINTERACTIONVOLUMESEEDS_HH
#define DUMUX_DISCRETIZATION_MPFA_O_HYBRIDFPS_GLOBALINTERACTIONVOLUMESEEDS_HH

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
class CCMpfaOHybridFpsGlobalInteractionVolumeSeeds
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Helper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolumeOMethod = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using OSeed = typename InteractionVolumeOMethod::Seed;
    using InteractionVolumeFpsMethod = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using FpsSeed = typename InteractionVolumeFpsMethod::Seed;
    using Element = typename GridView::template Codim<0>::Entity;

    using IndexType = typename GridView::IndexSet::IndexType;
    using LocalIndexType = typename InteractionVolumeOMethod::LocalIndexType;

    static const int dim = GridView::dimension;
    static const int numEq = GET_PROP_VALUE(TypeTag, NumEq);

public:
    CCMpfaOHybridFpsGlobalInteractionVolumeSeeds(const GridView gridView) : gridView_(gridView) {}

    // initializes the interaction volumes or the seeds
    void update(const Problem& problem, const std::vector<bool>& fpsVertices, const std::vector<bool>& boundaryVertices)
    {
        problemPtr_ = &problem;
        fpsSeeds_.clear();
        oSeeds_.clear();

        // -1 indicates that the scvf has not been handled yet
        auto numScvf = problem_().model().globalFvGeometry().numScvf();
        scvfIndexMap_.resize(numScvf, -1);

        // detect and handle the boundary first
        initializeFpsSeeds_(fpsVertices, boundaryVertices);
        initializeOSeeds_();
    }

    const OSeed& seed(const SubControlVolumeFace& scvf) const
    { return oSeeds_[scvfIndexMap_[scvf.index()]]; }

    const FpsSeed& boundarySeed(const SubControlVolumeFace& scvf, const LocalIndexType eqIdx) const
    { return fpsSeeds_[scvfIndexMap_[scvf.index()]]; }

private:
    void initializeFpsSeeds_(const std::vector<bool>& fpsVertices, const std::vector<bool>& boundaryVertices)
    {
        fpsSeeds_.reserve(fpsVertices.size());

        IndexType fpsSeedIndex = 0;
        for (const auto& element : elements(gridView_))
        {
            auto fvGeometry = localView(problem_().model().globalFvGeometry());
            fvGeometry.bind(element);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                // skip the rest if we already handled this face or if face doesn't belong to an fps element
                if (scvfIndexMap_[scvf.index()] != -1 || !fpsVertices[scvf.vertexIndex()])
                    continue;

                // the fps interaction volume seed
                auto fpsSeed = Helper::makeBoundaryInteractionVolumeSeed(problem_(), element, fvGeometry, scvf, boundaryVertices[scvf.vertexIndex()]);

                // update the index map entries for the global scv faces in the interaction volume
                for (const auto& localScvfSeed : fpsSeed.scvfSeeds())
                    for (const auto scvfIdxGlobal : localScvfSeed.globalScvfIndices())
                        scvfIndexMap_[scvfIdxGlobal] = fpsSeedIndex;

                // store interaction volume and increment counter
                fpsSeeds_.emplace_back(std::move(fpsSeed));
                fpsSeedIndex++;
            }
        }

        // shrink fps seed vector to actual size
        fpsSeeds_.shrink_to_fit();
    }

    void initializeOSeeds_()
    {
        oSeeds_.reserve(problem_().gridView().size(dim));

        IndexType oSeedIndex = 0;
        for (const auto& element : elements(gridView_))
        {
            auto fvGeometry = localView(problem_().model().globalFvGeometry());
            fvGeometry.bind(element);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (scvfIndexMap_[scvf.index()] != -1)
                    continue;

                // the inner interaction volume seed
                auto seed = Helper::makeInnerInteractionVolumeSeed(problem_(), element, fvGeometry, scvf);

                // update the index map entries for the global scv faces in the interaction volume
                for (const auto& localScvf : seed.scvfSeeds())
                    for (const auto scvfIdxGlobal : localScvf.globalScvfIndices())
                        scvfIndexMap_[scvfIdxGlobal] = oSeedIndex;

                // store interaction volume and increment counter
                oSeeds_.emplace_back(std::move(seed));
                oSeedIndex++;
            }
        }

        // shrink seed vector to actual size
        oSeeds_.shrink_to_fit();
    }

    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;
    GridView gridView_;
    std::vector<IndexType> scvfIndexMap_;
    std::vector<FpsSeed> fpsSeeds_;
    std::vector<OSeed> oSeeds_;
};
} // end namespace


#endif
