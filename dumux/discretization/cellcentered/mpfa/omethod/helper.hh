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
 * \brief Helper class to get the required information on an interaction volume.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFAO_HELPER_HH
#define DUMUX_DISCRETIZATION_CC_MPFAO_HELPER_HH

#include <dumux/common/math.hh>
#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

#include "localsubcontrolentities.hh"
#include "localsubcontrolentityseeds.hh"

namespace Dumux
{
/*!
 * \ingroup Mpfa
 * \brief Helper class to get the required information on an interaction volume.
 *        Specialization for the Mpfa-O method in two dimensions.
 */
template<class TypeTag>
class MpfaHelperBase<TypeTag, MpfaMethods::oMethod, 2>
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using LocalIndexType = typename InteractionVolume::LocalIndexType;

    using InteractionVolumeSeed = typename InteractionVolume::Seed;
    using ScvSeed = typename InteractionVolumeSeed::LocalScvSeed;
    using ScvfSeed = typename InteractionVolumeSeed::LocalScvfSeed;

    using GlobalIndexSet = std::vector<IndexType>;
    using LocalIndexSet = std::vector<LocalIndexType>;

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
public:
    static InteractionVolumeSeed makeInnerInteractionVolumeSeed(const Problem& problem,
                                                                const Element& element,
                                                                const FVElementGeometry& fvGeometry,
                                                                const SubControlVolumeFace& scvf)
    {
        std::vector<ScvSeed> scvSeeds;
        std::vector<ScvfSeed> scvfSeeds;

        fillEntitySeeds_(scvSeeds, scvfSeeds, problem, element, fvGeometry, scvf, /*dummy*/0);
        return InteractionVolumeSeed(std::move(scvSeeds), std::move(scvfSeeds), false);
    }

    static InteractionVolumeSeed makeBoundaryInteractionVolumeSeed(const Problem& problem,
                                                                   const Element& element,
                                                                   const FVElementGeometry& fvGeometry,
                                                                   const SubControlVolumeFace& scvf,
                                                                   const LocalIndexType eqIdx)
    {
        std::vector<ScvSeed> scvSeeds;
        std::vector<ScvfSeed> scvfSeeds;

        fillEntitySeeds_(scvSeeds, scvfSeeds, problem, element, fvGeometry, scvf, eqIdx);
        return InteractionVolumeSeed(std::move(scvSeeds), std::move(scvfSeeds), true);
    }

private:
    static void fillEntitySeeds_(std::vector<ScvSeed>& scvSeeds,
                                 std::vector<ScvfSeed>& scvfSeeds,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const SubControlVolumeFace& scvf,
                                 const LocalIndexType eqIdx)
    {
        // The vertex index around which we construct the interaction volume
        auto vIdxGlobal = scvf.vertexIndex();
        bool onBoundary = problem.model().globalFvGeometry().isBoundaryVertex(vIdxGlobal);

        // Get the two scv faces in the first scv
        auto scvfVector = Implementation::getScvFacesAtVertex(vIdxGlobal, element, fvGeometry);
        GlobalIndexSet globalScvfIndices({scvfVector[0]->index(), scvfVector[1]->index()});

        // The global index of the first scv of the interaction region
        auto scvIdx0 = scvf.insideScvIdx();

        // rotate counter clockwise and create the entities
        performRotation_(problem, scvfVector, scvSeeds, scvfSeeds, scvIdx0, eqIdx);
        if (onBoundary)
        {
            // the local scvf index of the second local face of the first local scv
            LocalIndexType storeIdx = scvfSeeds.size();

            // clockwise rotation until hitting the boundary again
            performRotation_(problem, scvfVector, scvSeeds, scvfSeeds, scvIdx0, eqIdx, /*clockwise*/true);

            // Finish by creating the first scv
            scvSeeds.emplace(scvSeeds.begin(), ScvSeed(std::move(globalScvfIndices),
                                                       LocalIndexSet({0, storeIdx}),
                                                       scvIdx0));
        }
        else
            // Finish by creating the first scv
            scvSeeds.emplace(scvSeeds.begin(), ScvSeed(std::move(globalScvfIndices),
                                                       LocalIndexSet({0, static_cast<LocalIndexType>(scvfSeeds.size()-1)}),
                                                       scvIdx0));
    }

    // clockwise rotation and construction of local scv & scv face entities
    template<class ScvfPointerVector>
    static void performRotation_(const Problem& problem,
                                 const ScvfPointerVector& scvfVector,
                                 std::vector<ScvSeed>& scvSeeds,
                                 std::vector<ScvfSeed>& scvfSeeds,
                                 const IndexType scvIdx0,
                                 const LocalIndexType eqIdx,
                                 const bool clockWise = false)
    {
        // extract the actual local indices from the containers
        IndexType localScvIdx = scvSeeds.size();
        IndexType localScvfIdx = scvfSeeds.size();

        // fvGeometry object to bind the neighbouring element during rotation
        auto outsideFvGeometry = localView(problem.model().globalFvGeometry());

        // Start/continue interaction region construction from the given scv face
        IndexType startScvfIdx = clockWise ? 1 : 0;
        auto curScvf = *scvfVector[startScvfIdx];
        bool firstIteration = true;
        bool finished = false;

        while (!finished)
        {
            // Get some indices beforehand
            IndexType globalScvfIdx = curScvf.index();
            IndexType insideGlobalScvIdx = curScvf.insideScvIdx();
            IndexType outsideGlobalScvIdx = curScvf.outsideScvIdx();
            IndexType insideLocalScvIdx = firstIteration ? 0 : localScvIdx;

            // the current element inside of the scv face
            auto insideElement = problem.model().globalFvGeometry().element(insideGlobalScvIdx);
            auto faceType = Implementation::getMpfaFaceType(problem, insideElement, curScvf, eqIdx);

            // if the face touches the boundary, create a boundary scvf entity
            if (curScvf.boundary())
            {
                assert(faceType == MpfaFaceTypes::neumann || faceType == MpfaFaceTypes::dirichlet);
                scvfSeeds.emplace_back(ScvfSeed(curScvf,
                                                LocalIndexSet({insideLocalScvIdx}),
                                                GlobalIndexSet({globalScvfIdx}),
                                                faceType));
                // rotation loop is finished
                finished = true; return;
            }

            // if outside scv is the first one again, finish loop
            if (outsideGlobalScvIdx == scvIdx0)
            {
                // create scv face entity for the last face of the loop
                scvfSeeds.emplace_back(ScvfSeed(curScvf,
                                                LocalIndexSet({insideLocalScvIdx, 0}),
                                                GlobalIndexSet({globalScvfIdx, scvfVector[1]->index()}),
                                                faceType));

                // create duplicate face if it is an interior boundary,
                if (faceType != MpfaFaceTypes::interior)
                    scvfSeeds.emplace_back(ScvfSeed(curScvf,
                                                    LocalIndexSet({0, insideLocalScvIdx}),
                                                    GlobalIndexSet({scvfVector[1]->index(), globalScvfIdx}),
                                                    faceType));

                // rotation loop is finished
                finished = true; return;
            }

            // If we get here, there are outside entities
            auto outsideElement = problem.model().globalFvGeometry().element(outsideGlobalScvIdx);
            outsideFvGeometry.bind(outsideElement);

            // get the two scv faces in the outside element that share the vertex
            auto outsideScvfVector = Implementation::getCommonAndNextScvFace(curScvf, outsideFvGeometry, clockWise);
            IndexType commonFaceCoordIdx = clockWise ? 0 : 1;
            IndexType nextFaceCoordIdx = clockWise ? 1 : 0;
            auto& commonScvf = *outsideScvfVector[commonFaceCoordIdx];
            auto& nextScvf = *outsideScvfVector[nextFaceCoordIdx];

            // create local scv face entity of the current scvf
            IndexType commonGlobalScvfIdx = commonScvf.index();
            IndexType outsideLocalScvIdx = outsideGlobalScvIdx == scvIdx0 ? 0 : localScvIdx+1;

            scvfSeeds.emplace_back(ScvfSeed(curScvf,
                                            LocalIndexSet({insideLocalScvIdx, outsideLocalScvIdx}),
                                            GlobalIndexSet({globalScvfIdx, commonGlobalScvfIdx}),
                                            faceType));
            localScvfIdx++;

            // create duplicate face if it is an interior boundary,
            if (faceType != MpfaFaceTypes::interior)
            {
                // face from "outside" to "inside", local index set changes
                scvfSeeds.emplace_back(ScvfSeed(commonScvf,
                                                LocalIndexSet({outsideLocalScvIdx, insideLocalScvIdx}),
                                                GlobalIndexSet({commonGlobalScvfIdx, globalScvfIdx}),
                                                faceType));
                localScvfIdx++;
            }

            // create index set storing the two local scvf indices
            LocalIndexSet localScvfs(2);
            localScvfs[commonFaceCoordIdx] = localScvfIdx-1;
            localScvfs[nextFaceCoordIdx] = localScvfIdx;

            // create "outside" scv
            GlobalIndexSet globalScvfIndices({outsideScvfVector[0]->index(), outsideScvfVector[1]->index()});
            scvSeeds.emplace_back(ScvSeed(std::move(globalScvfIndices), std::move(localScvfs), curScvf.outsideScvIdx()));
            localScvIdx++;

            // create the next scvf in the following iteration
            curScvf = nextScvf;
            firstIteration = false;
        }
    }
};

} // end namespace

#endif
