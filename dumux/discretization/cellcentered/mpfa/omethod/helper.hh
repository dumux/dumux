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
 *        Specialization for the Mpfa-O method in two dimensions embedded in a 2d world.
 */
template<class TypeTag>
class MpfaMethodHelper<TypeTag, MpfaMethods::oMethod, /*dim*/2, /*dimWorld*/2>
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, MpfaHelper);

    static const int dim = 2;
    static const int dimWorld = 2;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    // The o method is used on the boundaries also for other methods. These will use this helper
    // class to set up the seeds on the boundary. Therefore we extract all seed information from
    // the boundary interaction volume here to be compatible with other mpfa methods where seed types
    // most probably differ. We use the fact here, that for the o-method, boundary and interior
    // interaction volumes are identical, which is always given.
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using InteractionVolumeSeed = typename InteractionVolume::Seed;
    using ScvSeed = typename InteractionVolumeSeed::LocalScvSeed;
    using ScvfSeed = typename InteractionVolumeSeed::LocalScvfSeed;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalIndexType = typename InteractionVolume::GlobalIndexType;
    using LocalIndexType = typename InteractionVolume::LocalIndexType;

    using GlobalIndexSet = typename InteractionVolume::GlobalIndexSet;
    using LocalIndexSet = typename InteractionVolume::LocalIndexSet;
public:
    static InteractionVolumeSeed makeInnerInteractionVolumeSeed(const Problem& problem,
                                                                const Element& element,
                                                                const FVElementGeometry& fvGeometry,
                                                                const SubControlVolumeFace& scvf)
    {
        std::vector<ScvSeed> scvSeeds;
        std::vector<ScvfSeed> scvfSeeds;

        // reserve sufficient memory
        scvSeeds.reserve(20);
        scvfSeeds.reserve(20);

        fillEntitySeeds_(scvSeeds, scvfSeeds, problem, element, fvGeometry, scvf);

        // shrink containers to necessary size
        scvSeeds.shrink_to_fit();
        scvfSeeds.shrink_to_fit();

        return InteractionVolumeSeed(std::move(scvSeeds), std::move(scvfSeeds), false);
    }

    static InteractionVolumeSeed makeBoundaryInteractionVolumeSeed(const Problem& problem,
                                                                   const Element& element,
                                                                   const FVElementGeometry& fvGeometry,
                                                                   const SubControlVolumeFace& scvf)
    {
        std::vector<ScvSeed> scvSeeds;
        std::vector<ScvfSeed> scvfSeeds;

        // reserve sufficient memory
        scvSeeds.reserve(20);
        scvfSeeds.reserve(20);

        fillEntitySeeds_(scvSeeds, scvfSeeds, problem, element, fvGeometry, scvf);

        // shrink containers to necessary size
        scvSeeds.shrink_to_fit();
        scvfSeeds.shrink_to_fit();

        return InteractionVolumeSeed(std::move(scvSeeds), std::move(scvfSeeds), true);
    }

private:
    template<typename ScvSeedType, typename ScvfSeedType>
    static void fillEntitySeeds_(std::vector<ScvSeedType>& scvSeeds,
                                 std::vector<ScvfSeedType>& scvfSeeds,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const SubControlVolumeFace& scvf)
    {
        // Check whether or not we are touching the boundary here
        bool onBoundary = problem.model().globalFvGeometry().scvfTouchesBoundary(scvf);

        // Get the two scv faces in the first scv
        auto scvfVector = Implementation::getScvFacesAtVertex(scvf.vertexIndex(), element, fvGeometry);

        // The global index of the first scv of the interaction region
        auto scvIdx0 = scvf.insideScvIdx();

        // rotate counter clockwise and create the entities
        performRotation_(problem, scvfVector, scvSeeds, scvfSeeds, scvIdx0);

        if (onBoundary)
        {
            // the local scvf index of the second local face of the first local scv
            LocalIndexType storeIdx = scvfSeeds.size();

            // clockwise rotation until hitting the boundary again
            performRotation_(problem, scvfVector, scvSeeds, scvfSeeds, scvIdx0, /*clockwise*/true);

            // Finish by creating the first scv
            scvSeeds.emplace(scvSeeds.begin(), GlobalIndexSet({scvfVector[0]->index(), scvfVector[1]->index()}),
                                               LocalIndexSet({0, storeIdx}),
                                               scvIdx0);
        }
        else
            // Finish by creating the first scv
            scvSeeds.emplace(scvSeeds.begin(), GlobalIndexSet({scvfVector[0]->index(), scvfVector[1]->index()}),
                                               LocalIndexSet({0, static_cast<LocalIndexType>(scvfSeeds.size()-1)}),
                                               scvIdx0);
    }

    // in 2d we can make use of knowing the basis orientation (right hand system)
    // that way we can find the interaction regions more rapidly by clockwise and
    // counter clockwise rotation and construction of local scv & scv face entities
    template<typename ScvSeedType, typename ScvfSeedType, typename ScvfPointerVector>
    static void performRotation_(const Problem& problem,
                                 const ScvfPointerVector& scvfVector,
                                 std::vector<ScvSeedType>& scvSeeds,
                                 std::vector<ScvfSeedType>& scvfSeeds,
                                 const GlobalIndexType scvIdx0,
                                 const bool clockWise = false)
    {
        // extract the actual local indices from the containers
        LocalIndexType localScvIdx = scvSeeds.size();
        LocalIndexType localScvfIdx = scvfSeeds.size();

        // fvGeometry object to bind the neighbouring element during rotation
        auto outsideFvGeometry = localView(problem.model().globalFvGeometry());

        // Start/continue interaction region construction from the given scv face
        LocalIndexType startScvfIdx = clockWise ? 1 : 0;
        auto curScvf = *scvfVector[startScvfIdx];
        bool firstIteration = true;
        bool finished = false;

        while (!finished)
        {
            // Get some indices beforehand
            GlobalIndexType outsideGlobalScvIdx = curScvf.outsideScvIdx();
            LocalIndexType insideLocalScvIdx = firstIteration ? 0 : localScvIdx;

            // the current element inside of the scv face
            auto insideElement = problem.model().globalFvGeometry().element(curScvf.insideScvIdx());
            auto faceType = Implementation::getMpfaFaceType(problem, insideElement, curScvf);

            // if the face touches the boundary, create a boundary scvf entity
            if (curScvf.boundary())
            {
                assert(faceType == MpfaFaceTypes::neumann || faceType == MpfaFaceTypes::dirichlet);
                scvfSeeds.emplace_back( curScvf,
                                        insideLocalScvIdx,
                                        LocalIndexSet(),
                                        GlobalIndexSet(),
                                        faceType );
                // rotation loop is finished
                finished = true; return;
            }

            // if outside scv is the first one again, finish loop
            if (outsideGlobalScvIdx == scvIdx0)
            {
                // create scv face entity for the last face of the loop
                scvfSeeds.emplace_back( curScvf,
                                        insideLocalScvIdx,
                                        LocalIndexSet({0}),
                                        GlobalIndexSet({scvfVector[1]->index()}),
                                        faceType );

                // rotation loop is finished
                finished = true; return;
            }

            // If we get here, there are outside entities
            auto outsideElement = problem.model().globalFvGeometry().element(outsideGlobalScvIdx);
            outsideFvGeometry.bindElement(outsideElement);

            // get the two scv faces in the outside element that share the vertex
            auto outsideScvfVector = Implementation::getCommonAndNextScvFace(curScvf, outsideFvGeometry, clockWise);
            GlobalIndexType commonFaceCoordIdx = clockWise ? 0 : 1;
            GlobalIndexType nextFaceCoordIdx = clockWise ? 1 : 0;
            auto& commonScvf = *outsideScvfVector[commonFaceCoordIdx];
            auto& nextScvf = *outsideScvfVector[nextFaceCoordIdx];

            // create local scv face entity of the current scvf
            LocalIndexType outsideLocalScvIdx = localScvIdx+1;
            scvfSeeds.emplace_back( curScvf,
                                    insideLocalScvIdx,
                                    LocalIndexSet({outsideLocalScvIdx}),
                                    GlobalIndexSet({commonScvf.index()}),
                                    faceType );
            localScvfIdx++;

            // create index set storing the two local scvf indices
            LocalIndexSet localScvfs(2);
            localScvfs[commonFaceCoordIdx] = localScvfIdx-1;
            localScvfs[nextFaceCoordIdx] = localScvfIdx;

            // create "outside" scv
            GlobalIndexSet globalScvfIndices({outsideScvfVector[0]->index(), outsideScvfVector[1]->index()});
            scvSeeds.emplace_back(std::move(globalScvfIndices), std::move(localScvfs), outsideGlobalScvIdx);
            localScvIdx++;

            // create the next scvf in the following iteration
            curScvf = nextScvf;
            firstIteration = false;
        }
    }
};

/*!
 * \ingroup Mpfa
 * \brief Helper class to get the required information on an interaction volume.
 *        Specialization for the Mpfa-O method in two dimensions embedded in a 3d world.
 */
template<class TypeTag>
class MpfaMethodHelper<TypeTag, MpfaMethods::oMethod, /*dim*/2, /*dimWorld*/3>
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, MpfaHelper);

    static const int dim = 2;
    static const int dimWorld = 3;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    // The o method is used on the boundaries also for other methods. These will use this helper
    // class to set up the seeds on the boundary. Therefore we extract all seed information from
    // the boundary interaction volume here to be compatible with other mpfa methods where seed types
    // most probably differ. We use the fact here, that for the o-method, boundary and interior
    // interaction volumes are identical, which is always given.
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using InteractionVolumeSeed = typename InteractionVolume::Seed;
    using ScvSeed = typename InteractionVolumeSeed::LocalScvSeed;
    using ScvfSeed = typename InteractionVolumeSeed::LocalScvfSeed;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalIndexType = typename InteractionVolume::GlobalIndexType;
    using LocalIndexType = typename InteractionVolume::LocalIndexType;

    using GlobalIndexSet = typename InteractionVolume::GlobalIndexSet;
    using LocalIndexSet = typename InteractionVolume::LocalIndexSet;
public:
    static InteractionVolumeSeed makeInnerInteractionVolumeSeed(const Problem& problem,
                                                                const Element& element,
                                                                const FVElementGeometry& fvGeometry,
                                                                const SubControlVolumeFace& scvf)
    {
        // if the scvf does not touch a branching point, use simplified algorithm to create interaction volume seed
        if (!problem.model().globalFvGeometry().scvfTouchesBranchingPoint(scvf))
            return MpfaMethodHelper<TypeTag, MpfaMethods::oMethod, 2, 2>::makeInnerInteractionVolumeSeed(problem, element, fvGeometry, scvf);

        std::vector<ScvSeed> scvSeeds;
        std::vector<ScvfSeed> scvfSeeds;

        // reserve sufficient memory
        scvSeeds.reserve(20);
        scvfSeeds.reserve(20);

        // The vertex index around which we construct the interaction volume
        auto vIdxGlobal = scvf.vertexIndex();

        // Get the two scv faces in the scv
        auto scvfVector = Implementation::getScvFacesAtVertex(vIdxGlobal, element, fvGeometry);

        // fill the entity seed data
        fillEntitySeeds_(scvSeeds, scvfSeeds, problem, element, fvGeometry, scvfVector, vIdxGlobal);

        // shrink containers to necessary size
        scvSeeds.shrink_to_fit();
        scvfSeeds.shrink_to_fit();

        for (auto& scvfSeed : scvfSeeds)
            scvfSeed.makeOutsideDataUnique();

        return InteractionVolumeSeed(std::move(scvSeeds), std::move(scvfSeeds), false);
    }

    static InteractionVolumeSeed makeBoundaryInteractionVolumeSeed(const Problem& problem,
                                                                   const Element& element,
                                                                   const FVElementGeometry& fvGeometry,
                                                                   const SubControlVolumeFace& scvf)
    {
        // if the scvf does not touch a branching point, use simplified algorithm to create interaction volume seed
        if (!problem.model().globalFvGeometry().scvfTouchesBranchingPoint(scvf))
            return MpfaMethodHelper<TypeTag, MpfaMethods::oMethod, 2, 2>::makeBoundaryInteractionVolumeSeed(problem, element, fvGeometry, scvf);

        std::vector<ScvSeed> scvSeeds;
        std::vector<ScvfSeed> scvfSeeds;

        // reserve sufficient memory
        scvSeeds.reserve(20);
        scvfSeeds.reserve(20);

        // The vertex index around which we construct the interaction volume
        auto vIdxGlobal = scvf.vertexIndex();

        // Get the two scv faces in the scv
        auto scvfVector = Implementation::getScvFacesAtVertex(vIdxGlobal, element, fvGeometry);

        // fill the entity seed data
        fillEntitySeeds_(scvSeeds, scvfSeeds, problem, element, fvGeometry, scvfVector, vIdxGlobal);

        // shrink containers to necessary size
        scvSeeds.shrink_to_fit();
        scvfSeeds.shrink_to_fit();

        for (auto& scvfSeed : scvfSeeds)
            scvfSeed.makeOutsideDataUnique();

        return InteractionVolumeSeed(std::move(scvSeeds), std::move(scvfSeeds), true);
    }

private:
    template<class ScvSeedType, class ScvfSeedType, class ScvfVector>
    static void fillEntitySeeds_(std::vector<ScvSeedType>& scvSeeds,
                                 std::vector<ScvfSeedType>& scvfSeeds,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ScvfVector& scvfVector,
                                 const GlobalIndexType vIdxGlobal)
    {
        // make the scv without knowing the local scvf indices yet
        // take the inside scv index from the first scvf (is the same for all)
        scvSeeds.emplace_back( GlobalIndexSet({scvfVector[0]->index(),
                                               scvfVector[1]->index()}),
                               scvfVector[0]->insideScvIdx());

        // make the scvf seeds for the two scvfs connected to the scv
        auto& actualScvSeed = scvSeeds.back();
        LocalIndexType actualLocalScvIdx = scvSeeds.size()-1;

        for (int coordDir = 0; coordDir < dim; ++coordDir)
        {
            const auto& actualScvf = *scvfVector[coordDir];

            // if scvf is on a boundary, we create the scvfSeed and make no neighbor
            if (actualScvf.boundary())
            {
                // set the local scvfIndex of the face that is about to created
                actualScvSeed.setLocalScvfIndex(coordDir, scvfSeeds.size());

                // create the scvf seed
                scvfSeeds.emplace_back( actualScvf,
                                        actualLocalScvIdx,
                                        LocalIndexSet(),
                                        GlobalIndexSet(),
                                        Implementation::getMpfaFaceType(problem, element, actualScvf) );
            }
            else
            {
                // we loop over all neighbors of this face
                for (auto outsideGlobalScvIdx : actualScvf.outsideScvIndices())
                {
                    // get outside element, fvgeometry etc.
                    auto outsideElement = problem.model().globalFvGeometry().element(outsideGlobalScvIdx);
                    auto outsideFvGeometry = localView(problem.model().globalFvGeometry());
                    outsideFvGeometry.bindElement(outsideElement);

                    // find scvf in outside corresponding to the actual scvf
                    auto outsideScvfVector = Implementation::getScvFacesAtVertex(vIdxGlobal, outsideElement, outsideFvGeometry);
                    auto commonFaceLocalIdx = Implementation::getCommonFaceLocalIndex(actualScvf, outsideScvfVector);
                    const auto& outsideScvf = *outsideScvfVector[commonFaceLocalIdx];
                    const auto outsideScvfIdx = outsideScvf.index();

                    // check if the outside scv already exists and get its local index
                    bool outsideScvExists = false;
                    LocalIndexType outsideLocalScvIdx = 0;
                    for (auto&& scvSeed : scvSeeds)
                    {
                        if (scvSeed.globalIndex() == outsideGlobalScvIdx)
                        {
                            outsideScvExists = true;
                            break;
                        }
                        // keep track of local index
                        outsideLocalScvIdx++;
                    }

                    // check if face already exists inside or in any neighbor
                    bool outsideScvfExists = false;
                    bool insideScvfExists = false;
                    LocalIndexType outsideLocalScvfIdx = 0;
                    LocalIndexType insideLocalScvfIdx = 0;
                    for (auto&& scvfSeed : scvfSeeds)
                    {
                        if (Implementation::contains(actualScvf.outsideScvIndices(), scvfSeed.insideGlobalScvIndex()) &&
                            Implementation::contains(scvfSeed.outsideGlobalScvIndices(), actualScvf.insideScvIdx()))
                            outsideScvfExists = true;
                        else
                            if (!outsideScvfExists)
                                outsideLocalScvfIdx++;

                        if (scvfSeed.insideGlobalScvIndex() == actualScvf.insideScvIdx() &&
                            scvfSeed.outsideGlobalScvIndices().size() == actualScvf.outsideScvIndices().size() &&
                            std::equal(scvfSeed.outsideGlobalScvIndices().begin(), scvfSeed.outsideGlobalScvIndices().end(), actualScvf.outsideScvIndices().begin()))
                            insideScvfExists = true;
                        else
                            if (!insideScvfExists)
                                insideLocalScvfIdx++;
                    }

                    // we should never have two local scv faces existing inside and in an outside element
                    assert(!(outsideScvfExists && insideScvfExists) && "The scv face seems to exist twice!");

                    // outside scv has to be created
                    if (!outsideScvExists)
                    {
                        // The face does not yet exist
                        if (!insideScvfExists && !outsideScvfExists)
                        {
                            // set the local scvfIndex of the face that is about to created
                            actualScvSeed.setLocalScvfIndex(coordDir, scvfSeeds.size());

                            // create scvf seed
                            scvfSeeds.emplace_back(actualScvf,
                                                   actualLocalScvIdx,
                                                   Implementation::getMpfaFaceType(problem, element, actualScvf));

                            // pass the actual outside indices to the new scvf seed
                            scvfSeeds.back().addOutsideData(outsideScvfIdx, static_cast<LocalIndexType>(scvSeeds.size()));
                        }
                        else if (insideScvfExists && !outsideScvfExists)
                        {
                            // pass info on outside to the inside scvf seed
                            scvfSeeds[insideLocalScvfIdx].addOutsideData(outsideScvfIdx, static_cast<LocalIndexType>(scvSeeds.size()));
                        }
                        else if (!insideScvfExists && outsideScvfExists)
                        {
                            // set the local scvfIndex of the outside face
                            actualScvSeed.setLocalScvfIndex(coordDir, outsideLocalScvfIdx);

                            // pass info on inside to the outside scvf seed
                            scvfSeeds[outsideLocalScvfIdx].addOutsideData(actualScvf.index(), actualLocalScvIdx);
                        }

                        // make outside scv by recursion
                        fillEntitySeeds_(scvSeeds, scvfSeeds, problem, outsideElement, outsideFvGeometry, outsideScvfVector, vIdxGlobal);
                    }
                    else
                    {
                        // outside scv exists, but no corresponding local scv face yet. Make it from inside.
                        if (!outsideScvfExists && !insideScvfExists)
                        {
                            // set the local scvfIndex of the face that is about to created
                            actualScvSeed.setLocalScvfIndex(coordDir, scvfSeeds.size());

                            // create scvf seed
                            scvfSeeds.emplace_back(actualScvf,
                                                   actualLocalScvIdx,
                                                   Implementation::getMpfaFaceType(problem, element, actualScvf));

                            // pass the actual outside indices to the new scvf seed
                            scvfSeeds.back().addOutsideData(outsideScvfIdx, outsideLocalScvIdx);
                        }
                        else if (!insideScvfExists && outsideScvfExists)
                        {
                            // set the local scvfIndex of the found outside local scv face seed
                            actualScvSeed.setLocalScvfIndex(coordDir, outsideLocalScvfIdx);

                            // pass info on inside to the outside found local scvf seed
                            scvfSeeds[outsideLocalScvfIdx].addOutsideData(actualScvf.index(), actualLocalScvIdx);
                        }
                    }
                }
            }
        }
    }
};

/*!
 * \ingroup Mpfa
 * \brief Helper class to get the required information on an interaction volume.
 *        Specialization for the Mpfa-O method in three dimensions.
 */
template<class TypeTag>
class MpfaMethodHelper<TypeTag, MpfaMethods::oMethod, /*dim*/3, /*dimWorld*/3>
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, MpfaHelper);

    static const int dim = 3;
    static const int dimWorld = 3;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    // The o method is used on the boundaries also for other methods. These will use this helper
    // class to set up the seeds on the boundary. Therefore we extract all seed information from
    // the boundary interaction volume here to be compatible with other mpfa methods where seed types
    // most probably differ. We use the fact here, that for the o-method, boundary and interior
    // interaction volumes are identical, which is always given.
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using InteractionVolumeSeed = typename InteractionVolume::Seed;
    using ScvSeed = typename InteractionVolumeSeed::LocalScvSeed;
    using ScvfSeed = typename InteractionVolumeSeed::LocalScvfSeed;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalIndexType = typename InteractionVolume::GlobalIndexType;
    using LocalIndexType = typename InteractionVolume::LocalIndexType;

    using GlobalIndexSet = typename InteractionVolume::GlobalIndexSet;
    using LocalIndexSet = typename InteractionVolume::LocalIndexSet;

public:
    static InteractionVolumeSeed makeInnerInteractionVolumeSeed(const Problem& problem,
                                                                const Element& element,
                                                                const FVElementGeometry& fvGeometry,
                                                                const SubControlVolumeFace& scvf)
    {
        std::vector<ScvSeed> scvSeeds;
        std::vector<ScvfSeed> scvfSeeds;

        // reserve sufficient memory
        scvSeeds.reserve(50);
        scvfSeeds.reserve(50);

        // The vertex index around which we construct the interaction volume
        auto vIdxGlobal = scvf.vertexIndex();

        // Get the three scv faces in the scv
        auto scvfVector = Implementation::getScvFacesAtVertex(vIdxGlobal, element, fvGeometry);

        // create the scv entity seeds
        fillEntitySeeds_(scvSeeds, scvfSeeds, problem, element, fvGeometry, scvfVector, vIdxGlobal);

        // shrink containers to necessary size
        scvSeeds.shrink_to_fit();
        scvfSeeds.shrink_to_fit();

        return InteractionVolumeSeed(std::move(scvSeeds), std::move(scvfSeeds), false);
    }

    static InteractionVolumeSeed makeBoundaryInteractionVolumeSeed(const Problem& problem,
                                                                   const Element& element,
                                                                   const FVElementGeometry& fvGeometry,
                                                                   const SubControlVolumeFace& scvf)
    {
        std::vector<ScvSeed> scvSeeds;
        std::vector<ScvfSeed> scvfSeeds;

        // reserve sufficient memory
        scvSeeds.reserve(50);
        scvfSeeds.reserve(50);

        // The vertex index around which we construct the interaction volume
        auto vIdxGlobal = scvf.vertexIndex();

        // Get the three scv faces in the scv
        auto scvfVector = Implementation::getScvFacesAtVertex(vIdxGlobal, element, fvGeometry);

        // create the scv entity seeds
        fillEntitySeeds_(scvSeeds, scvfSeeds, problem, element, fvGeometry, scvfVector, vIdxGlobal);

        // shrink containers to necessary size
        scvSeeds.shrink_to_fit();
        scvfSeeds.shrink_to_fit();

        return InteractionVolumeSeed(std::move(scvSeeds), std::move(scvfSeeds), true);
    }

private:
    template<class ScvSeedType, class ScvfSeedType, class ScvfVector>
    static void fillEntitySeeds_(std::vector<ScvSeedType>& scvSeeds,
                                 std::vector<ScvfSeedType>& scvfSeeds,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ScvfVector& scvfVector,
                                 const GlobalIndexType vIdxGlobal)
    {
        // make the scv without knowing the local scvf indices yet
        // take the inside scv index from the first scvf (is the same for all)
        scvSeeds.emplace_back( GlobalIndexSet({scvfVector[0]->index(),
                                               scvfVector[1]->index(),
                                               scvfVector[2]->index()}),
                               scvfVector[0]->insideScvIdx());

        // make the scvf seeds for the three scvfs connected to the scv
        auto& actualScvSeed = scvSeeds.back();
        LocalIndexType actualLocalScvIdx = scvSeeds.size()-1;

        for (int coordDir = 0; coordDir < dim; ++coordDir)
        {
            const auto& actualScvf = *scvfVector[coordDir];

            // if scvf is on a boundary, we create the scvfSeed and make no neighbor
            if (actualScvf.boundary())
            {
                // set the local scvfIndex of the face that is about to created
                actualScvSeed.setLocalScvfIndex(coordDir, scvfSeeds.size());

                // create the scvf seed
                auto faceType = Implementation::getMpfaFaceType(problem, element, actualScvf);
                scvfSeeds.emplace_back( actualScvf,
                                        actualLocalScvIdx,
                                        LocalIndexSet(),
                                        GlobalIndexSet(),
                                        faceType );
            }
            else
            {
                auto outsideGlobalScvIdx = actualScvf.outsideScvIdx();
                auto globalScvfIndex = actualScvf.index();

                // check if the outside scv already exists and get its local index
                bool outsideExists = false;
                LocalIndexType outsideLocalScvIdx = 0;
                for (auto&& scvSeed : scvSeeds)
                {
                    if (scvSeed.globalIndex() == outsideGlobalScvIdx)
                    {
                        outsideExists = true;
                        break;
                    }
                    // keep track of local index
                    outsideLocalScvIdx++;
                }

                // if outside scv does not exist we have to make the scvf and the outside scv
                if (!outsideExists)
                {
                    // set the local scvfIndex of the face that is about to created
                    actualScvSeed.setLocalScvfIndex(coordDir, scvfSeeds.size());

                    // get outside element, fvgeometry etc.
                    auto outsideElement = problem.model().globalFvGeometry().element(outsideGlobalScvIdx);
                    auto outsideFvGeometry = localView(problem.model().globalFvGeometry());
                    outsideFvGeometry.bindElement(outsideElement);

                    // find scvf in outside corresponding to the actual scvf
                    auto outsideScvfVector = Implementation::getScvFacesAtVertex(vIdxGlobal, outsideElement, outsideFvGeometry);
                    auto commonFaceLocalIdx = Implementation::getCommonFaceLocalIndex(actualScvf, outsideScvfVector);
                    const auto& outsideScvf = *outsideScvfVector[commonFaceLocalIdx];

                    // create scvf seed
                    scvfSeeds.emplace_back(actualScvf,
                                           actualLocalScvIdx,
                                           LocalIndexSet( {static_cast<LocalIndexType>(scvSeeds.size())} ),
                                           GlobalIndexSet( {outsideScvf.index()} ),
                                           Implementation::getMpfaFaceType(problem, element, actualScvf));

                    // make outside scv by recursion
                    fillEntitySeeds_(scvSeeds, scvfSeeds, problem, outsideElement, outsideFvGeometry, outsideScvfVector, vIdxGlobal);
                }
                // we have to find out if it is necessary to create a new scvf
                else
                {
                    // find the scvf seed with the actual scvf as outside scvf
                    bool found = false;
                    LocalIndexType localScvfIdx = 0;
                    for (auto&& scvfSeed : scvfSeeds)
                    {
                        // boundary scvf seeds have no outside scvf
                        if (!scvfSeed.boundary() && scvfSeed.outsideGlobalScvfIndex() == globalScvfIndex)
                        {
                            // pass local scvf index to local scv
                            actualScvSeed.setLocalScvfIndex(coordDir, localScvfIdx);

                            // we found the corresponding face
                            found = true; break;
                        }
                        // keep track of local index
                        localScvfIdx++;
                    }

                    // if no corresponding scvf has been found, create it
                    if (!found)
                    {
                        // set the local scvfIndex of the face that is about to created
                        actualScvSeed.setLocalScvfIndex(coordDir, scvfSeeds.size());

                        // get outside element, fvgeometry etc.
                        auto outsideElement = problem.model().globalFvGeometry().element(outsideGlobalScvIdx);
                        auto outsideFvGeometry = localView(problem.model().globalFvGeometry());
                        outsideFvGeometry.bindElement(outsideElement);

                        // find scvf in outside corresponding to the actual scvf
                        auto outsideScvfVector = Implementation::getScvFacesAtVertex(vIdxGlobal, outsideElement, outsideFvGeometry);
                        auto commonFaceLocalIdx = Implementation::getCommonFaceLocalIndex(actualScvf, outsideScvfVector);
                        const auto& outsideScvf = *outsideScvfVector[commonFaceLocalIdx];

                        // make scv face seed
                        scvfSeeds.emplace_back(actualScvf,
                                               actualLocalScvIdx,
                                               LocalIndexSet( {outsideLocalScvIdx} ),
                                               GlobalIndexSet( {outsideScvf.index()} ),
                                               Implementation::getMpfaFaceType(problem, element, actualScvf));
                    }
                }
            }
        }
    }
};

} // end namespace

#endif
