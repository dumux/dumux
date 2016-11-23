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
#ifndef DUMUX_DISCRETIZATION_CC_MPFAO_HYBRID_FPS_HELPER_HH
#define DUMUX_DISCRETIZATION_CC_MPFAO_HYBRID_FPS_HELPER_HH

#include <dumux/discretization/cellcentered/mpfa/omethod/helper.hh>

namespace Dumux
{
/*!
 * \ingroup Mpfa
 * \brief Helper class to get the required information on an interaction volume.
 *        Specialization for the Mpfa-O method in two dimensions.
 */
template<class TypeTag, int dim>
class MpfaHelperBase<TypeTag, MpfaMethods::oMethodFpsHybrid, dim> : public MpfaHelperBase<TypeTag, MpfaMethods::oMethod, dim>
{
    using Base = MpfaHelperBase<TypeTag, MpfaMethods::oMethod, dim>;

    using Implementation = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);

    using GlobalIndexType = typename Base::GlobalIndexType;
    using LocalIndexType = typename Base::LocalIndexType;
    using GlobalIndexSet = typename Base::GlobalIndexSet;
    using LocalIndexSet = typename Base::LocalIndexSet;

    using InteractionVolumeSeed = typename InteractionVolume::Seed;
    using ScvSeed = typename InteractionVolumeSeed::LocalScvSeed;
    using ScvfSeed = typename InteractionVolumeSeed::LocalScvfSeed;

    using Element = typename GridView::template Codim<0>::Entity;

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
                                                                   const bool isBoundary)
    {
        std::vector<ScvSeed> scvSeeds;
        std::vector<ScvfSeed> scvfSeeds;

        fillEntitySeeds_(scvSeeds, scvfSeeds, problem, element, fvGeometry, scvf, isBoundary);
        return InteractionVolumeSeed(std::move(scvSeeds), std::move(scvfSeeds), isBoundary);
    }

private:
    static void fillEntitySeeds_(std::vector<ScvSeed>& scvSeeds,
                                 std::vector<ScvfSeed>& scvfSeeds,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const SubControlVolumeFace& scvf,
                                 const bool isBoundary)
    {
        // Get the two scv faces in the first scv
        auto scvfVector = Implementation::getScvFacesAtVertex(scvf.vertexIndex(), element, fvGeometry);

        // The global index of the first scv of the interaction region
        auto scvIdx0 = scvf.insideScvIdx();

        // rotate counter clockwise and create the entities
        Implementation::performRotation_(problem, scvfVector, scvSeeds, scvfSeeds, scvIdx0, /*dummy*/0);

        if (isBoundary)
        {
            // the local scvf index of the second local face of the first local scv
            LocalIndexType storeIdx = scvfSeeds.size();

            // clockwise rotation until hitting the boundary again
            Implementation::performRotation_(problem, scvfVector, scvSeeds, scvfSeeds, scvIdx0, /*dummy*/0, /*clockwise*/true);

            // Finish by creating the first scv
            scvSeeds.emplace(scvSeeds.begin(), ScvSeed(GlobalIndexSet({scvfVector[0]->index(), scvfVector[1]->index()}),
                                                       LocalIndexSet({0, storeIdx}),
                                                       scvIdx0));
        }
        else
            // Finish by creating the first scv
            scvSeeds.emplace(scvSeeds.begin(), ScvSeed(GlobalIndexSet({scvfVector[0]->index(), scvfVector[1]->index()}),
                                                       LocalIndexSet({0, static_cast<LocalIndexType>(scvfSeeds.size()-1)}),
                                                       scvIdx0));
    }
};

} // end namespace

#endif
