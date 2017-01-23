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
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_L_HELPER_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_L_HELPER_HH

#include <dumux/common/math.hh>
#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

#include "localsubcontrolentityseeds.hh"

namespace Dumux
{
/*!
 * \ingroup Mpfa
 * \brief Helper class to get the required information on an interaction volume.
 *        Specialization for the Mpfa-L method in two dimensions.
 */
template<class TypeTag>
class MpfaMethodHelper<TypeTag, MpfaMethods::lMethod, /*dim*/2, /*dimWorld*/2>
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, MpfaHelper);

    static const int dim = 2;
    static const int dimWorld = 2;

    // The mpfa-o helper class used to construct the boundary interaction volume seeds
    using oMethodHelper = CCMpfaHelperImplementation<TypeTag, MpfaMethods::oMethod, dim, dimWorld>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);

    using Element = typename GridView::template Codim<0>::Entity;

    using InteractionVolumeSeed = typename InteractionVolume::Seed;
    using ScvSeed = typename InteractionVolumeSeed::LocalScvSeed;
    using OuterScvSeed = typename InteractionVolumeSeed::LocalOuterScvSeed;
    using BoundaryInteractionVolumeSeed = typename BoundaryInteractionVolume::Seed;

    using GlobalIndexSet = typename InteractionVolume::GlobalIndexSet;
    using GlobalIndexType = typename InteractionVolume::GlobalIndexType;
    using LocalIndexType = typename InteractionVolume::LocalIndexType;
    using Matrix = typename InteractionVolume::Matrix;
public:
    static InteractionVolumeSeed makeInnerInteractionVolumeSeed(const Problem& problem,
                                                                const Element& element,
                                                                const FVElementGeometry& fvGeometry,
                                                                const SubControlVolumeFace& scvf)
    {
        std::vector<ScvSeed> scvSeeds;
        std::vector<OuterScvSeed> outerScvSeeds;
        std::vector<GlobalIndexType> globalScvfIndices(2);

        // we'll have maximal 2 ScvSeeds and 2 OuterScvSeeds
        scvSeeds.reserve(2);
        outerScvSeeds.reserve(2);

        fillEntitySeeds_(scvSeeds, outerScvSeeds, globalScvfIndices, problem, element, fvGeometry, scvf);

        // return interaction volume seed
        return InteractionVolumeSeed(std::move(scvSeeds), std::move(outerScvSeeds), std::move(globalScvfIndices));
    }

    static BoundaryInteractionVolumeSeed makeBoundaryInteractionVolumeSeed(const Problem& problem,
                                                                           const Element& element,
                                                                           const FVElementGeometry& fvGeometry,
                                                                           const SubControlVolumeFace& scvf)
    { return oMethodHelper::makeBoundaryInteractionVolumeSeed(problem, element, fvGeometry, scvf); }

    template<class InteractionRegion>
    static LocalIndexType selectionCriterion(const InteractionRegion& I1,
                                             const InteractionRegion& I2,
                                             const Matrix& M1,
                                             const Matrix& M2)
    {
        Scalar eps = 1e-10;
        Scalar t11, t12, t21, t22;
        t11 = M1[I1.contiFaceLocalIdx][0];
        t21 = M2[I2.contiFaceLocalIdx][0];
        if (I1.contiFaceLocalIdx == 0)
        {
            t12 = M1[I1.contiFaceLocalIdx][1];
            t22 = M2[I2.contiFaceLocalIdx][2];
        }
        else
        {
            t12 = M1[I1.contiFaceLocalIdx][2];
            t22 = M2[I2.contiFaceLocalIdx][1];
        }

        Scalar s1 = std::abs(t11-t12);
        Scalar s2 = std::abs(t22-t21);

        if (s1 < s2 + eps*s1)
            return 0;
        else
            return 1;
    }

private:
    static void fillEntitySeeds_(std::vector<ScvSeed>& scvSeeds,
                                 std::vector<OuterScvSeed>& outerScvSeeds,
                                 std::vector<GlobalIndexType>& globalScvfIndices,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const SubControlVolumeFace& scvf)
    {
        // make the first scv seed, we know this element will NOT be on the lowest local level
        auto scvfVector = Implementation::getScvFacesAtVertex(scvf.vertexIndex(), element, fvGeometry);
        auto localScvfIdx = Implementation::getLocalFaceIndex(scvf, scvfVector);
        scvSeeds.emplace_back( GlobalIndexSet({scvfVector[0]->index(), scvfVector[1]->index()}),
                               scvf.insideScvIdx(),
                               localScvfIdx );

        // get the surrounding elements and "outside" data
        LocalIndexType otherScvfIdx = 1-localScvfIdx;
        auto e2 = problem.model().globalFvGeometry().element(scvf.outsideScvIdx());
        auto e3 = problem.model().globalFvGeometry().element(scvfVector[otherScvfIdx]->outsideScvIdx());

        auto e2Geometry = localView(problem.model().globalFvGeometry());
        auto e3Geometry = localView(problem.model().globalFvGeometry());

        e2Geometry.bindElement(e2);
        e3Geometry.bindElement(e3);

        auto e2Scvfs = Implementation::getCommonAndNextScvFace(scvf, e2Geometry, /*clockwise?*/localScvfIdx == 1);
        auto e3Scvfs = Implementation::getCommonAndNextScvFace(*scvfVector[otherScvfIdx], e3Geometry, /*clockwise?*/localScvfIdx == 0);

        // we now know the two faces for which flux calculation will happen using this iv seed
        globalScvfIndices[0] = scvf.index();
        globalScvfIndices[1] = e2Scvfs[otherScvfIdx]->index();

        // scv seed for e2, we know the local common scvf index will be otherScvfIdx in 2d
        scvSeeds.emplace_back( GlobalIndexSet({e2Scvfs[0]->index(), e2Scvfs[1]->index()}),
                               e2Scvfs[0]->insideScvIdx(),
                               otherScvfIdx );

        // Outer seed for e3, we know the local common scvf index will be localScvfIdx in 2d
        outerScvSeeds.emplace_back(e3Scvfs[localScvfIdx]->insideScvIdx(),
                                   e3Scvfs[localScvfIdx]->index());

        // Outer seed for outside of e2, we know the local scvf index here will be localScvfIdx in 2d
        auto e4 = problem.model().globalFvGeometry().element(e2Scvfs[localScvfIdx]->outsideScvIdx());
        auto e4Geometry = localView(problem.model().globalFvGeometry());
        e4Geometry.bindElement(e4);
        auto e4Scvfs = Implementation::getCommonAndNextScvFace(*e2Scvfs[localScvfIdx], e4Geometry, /*clockwise?*/localScvfIdx == 1);

        outerScvSeeds.emplace_back(e4Scvfs[otherScvfIdx]->insideScvIdx(),
                                   e4Scvfs[otherScvfIdx]->index());
    }
};

/*!
 * \ingroup Mpfa
 * \brief Helper class to get the required information on an interaction volume.
 *        Specialization for the Mpfa-L method for 2d embedded in 3d space.
 */
template<class TypeTag>
class MpfaMethodHelper<TypeTag, MpfaMethods::lMethod, /*dim*/2, /*dimWorld*/3>
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, MpfaHelper);

    static const int dim = 2;
    static const int dimWorld = 3;

    // The mpfa-o helper class used to construct the boundary interaction volume seeds
    using oMethodHelper = CCMpfaHelperImplementation<TypeTag, MpfaMethods::oMethod, dim, dimWorld>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);

    using Element = typename GridView::template Codim<0>::Entity;

    using InteractionVolumeSeed = typename InteractionVolume::Seed;
    using BoundaryInteractionVolumeSeed = typename BoundaryInteractionVolume::Seed;

    using LocalIndexType = typename InteractionVolume::LocalIndexType;
    using Matrix = typename InteractionVolume::Matrix;
public:

    //! We know that this is only called for scvfs NOT touching:
    //!     - branching points
    //!     - domain boundaries
    //!     - interior boundaries
    //! Thus, here we can use the same algorithm as in the 2d case.
    static InteractionVolumeSeed makeInnerInteractionVolumeSeed(const Problem& problem,
                                                                const Element& element,
                                                                const FVElementGeometry& fvGeometry,
                                                                const SubControlVolumeFace& scvf)
    {
        return MpfaMethodHelper<TypeTag, MpfaMethods::lMethod, 2, 2>::makeInnerInteractionVolumeSeed(problem,
                                                                                                     element,
                                                                                                     fvGeometry,
                                                                                                     scvf);
    }

    //! Here we simply forward to the o-method helper (we use o-interaction volumes on the boundary)
    static BoundaryInteractionVolumeSeed makeBoundaryInteractionVolumeSeed(const Problem& problem,
                                                                           const Element& element,
                                                                           const FVElementGeometry& fvGeometry,
                                                                           const SubControlVolumeFace& scvf)
    { return oMethodHelper::makeBoundaryInteractionVolumeSeed(problem, element, fvGeometry, scvf); }

    //! The selection criterion is the same as in the dim = 2 = dimWorld case
    template<class InteractionRegion>
    static LocalIndexType selectionCriterion(const InteractionRegion& I1,
                                             const InteractionRegion& I2,
                                             const Matrix& M1,
                                             const Matrix& M2)
    { return MpfaMethodHelper<TypeTag, MpfaMethods::lMethod, 2, 2>::selectionCriterion(I1, I2, M1, M2); }
};

} // end namespace

#endif
