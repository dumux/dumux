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
 * \brief Base class for the global interaction volume seeds of mpfa methods.
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_L_GLOBALINTERACTIONVOLUMESEEDS_HH
#define DUMUX_DISCRETIZATION_MPFA_L_GLOBALINTERACTIONVOLUMESEEDS_HH

#include <dumux/discretization/cellcentered/mpfa/globalinteractionvolumeseedsbase.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

namespace Dumux
{
/*!
 * \ingroup Mpfa
 * \brief Specialization of the class for the mpfa-l method.
 */
template<class TypeTag>
class CCMpfaGlobalInteractionVolumeSeedsImplementation<TypeTag, MpfaMethods::lMethod>
       : public CCMpfaGlobalInteractionVolumeSeedsBase<TypeTag>
{
    using ParentType = CCMpfaGlobalInteractionVolumeSeedsBase<TypeTag>;

    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Helper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);

    using IndexType = typename InteractionVolume::GlobalIndexSet::value_type;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    CCMpfaGlobalInteractionVolumeSeedsImplementation(const GridView& gridView) : ParentType(gridView) {}

    template<typename SeedVector, typename BoundarySeedVector>
    void initializeSeeds(const std::vector<bool>& interiorOrDomainBoundaryVertices,
                         std::vector<IndexType>& scvfIndexMap,
                         SeedVector& seeds,
                         BoundarySeedVector& boundarySeeds)
    {
        seeds.clear();
        boundarySeeds.clear();
        scvfIndexMap.clear();

        // reserve memory
        const auto numScvf = this->problem().model().fvGridGeometry().numScvf();
        const auto numIBScvf = this->problem().model().fvGridGeometry().numInteriorBoundaryScvf();
        const auto numDBScvf = this->problem().model().fvGridGeometry().numDomainBoundaryScvf();
        const auto numBPScvf = this->problem().model().fvGridGeometry().numBranchingPointScvf();

        const int numLMethodIVs = (numScvf-numIBScvf-numDBScvf-numBPScvf)/2;
        const auto numOMethodIVs = this->problem().model().fvGridGeometry().numInteriorOrDomainBoundaryVertices()
                                   + this->problem().model().fvGridGeometry().numBranchingPointVertices();

        if (numLMethodIVs > 0)
            seeds.reserve( numLMethodIVs );
        boundarySeeds.reserve(numOMethodIVs);
        scvfIndexMap.resize(numScvf);

        // Keep track of which faces have been handled already
        std::vector<bool> isFaceHandled(numScvf, false);

        IndexType boundarySeedIndex = 0;
        IndexType seedIndex = 0;
        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(this->problem().model().fvGridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scvf : scvfs(fvGeometry))
            {
                // skip the rest if we already handled this face
                if (isFaceHandled[scvf.index()])
                    continue;

                // make boundary interaction volume seeds (o-method seeds) on:
                //      - interior boundaries
                //      - domain boundaries
                //      - branching points (l-method can not handle this)
                if (interiorOrDomainBoundaryVertices[scvf.vertexIndex()]
                    || this->problem().model().fvGridGeometry().touchesBranchingPoint(scvf))
                {
                    // make the boundary interaction volume seed
                    boundarySeeds.emplace_back(Helper::makeBoundaryInteractionVolumeSeed(this->problem(),
                                                                                         element,
                                                                                         fvGeometry,
                                                                                         scvf));

                    // update the index map entries for the global scv faces in the interaction volume
                    for (auto scvfIdxGlobal : boundarySeeds.back().globalScvfIndices())
                    {
                        scvfIndexMap[scvfIdxGlobal] = boundarySeedIndex;
                        isFaceHandled[scvfIdxGlobal] = true;
                    }

                    // increment counter
                    boundarySeedIndex++;
                }
                else
                {
                    // make the inner interaction volume seed only if we are on highest level of all connected elements
                    if (isLocalMaxLevel_(element, scvf))
                    {
                        seeds.emplace_back(Helper::makeInnerInteractionVolumeSeed(this->problem(),
                                                                                  element,
                                                                                  fvGeometry,
                                                                                  scvf));

                        // update the index map entries for the global scv faces in the interaction volume
                        for (auto scvfIdxGlobal : seeds.back().globalScvfIndices())
                        {
                            scvfIndexMap[scvfIdxGlobal] = seedIndex;
                            isFaceHandled[scvfIdxGlobal] = true;
                        }

                        // increment counter
                        seedIndex++;
                    }
                }
            }
        }
    }

    bool isLocalMaxLevel_(const Element& element, const SubControlVolumeFace& scvf) const
    { return this->problem().model().fvGridGeometry().element(scvf.outsideScvIdx()).level() <= element.level(); }
};

} // end namespace


#endif
