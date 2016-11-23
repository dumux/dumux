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
class CCMpfaGlobalInteractionVolumeSeedsImplementation<TypeTag, MpfaMethods::lMethod> : public CCMpfaGlobalInteractionVolumeSeedsBase<TypeTag>
{
    using ParentType = CCMpfaGlobalInteractionVolumeSeedsBase<TypeTag>;
    // the parent needs to be friend to access the private initialization method
    friend ParentType;

    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Helper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);

    using IndexType = typename InteractionVolume::GlobalIndexSet::value_type;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    CCMpfaGlobalInteractionVolumeSeedsImplementation(const GridView gridView) : ParentType(gridView)  {}

private:
    void initializeSeeds_(std::vector<bool>& boundaryVertices, std::vector<bool>& isFaceHandled)
    {
        const auto numScvf = this->problem_().model().globalFvGeometry().numScvf();
        const auto numBoundaryScvf = this->problem_().model().globalFvGeometry().numBoundaryScvf();

        // reserve memory
        this->seeds_.reserve(numScvf - numBoundaryScvf);
        this->boundarySeeds_.reserve(numBoundaryScvf);

        IndexType boundarySeedIndex = 0;
        IndexType seedIndex = 0;
        for (const auto& element : elements(this->gridView_))
        {
            auto fvGeometry = localView(this->problem_().model().globalFvGeometry());
            fvGeometry.bindElement(element);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                // skip the rest if we already handled this face
                if (isFaceHandled[scvf.index()])
                    continue;

                if (boundaryVertices[scvf.vertexIndex()])
                {
                    // make the boundary interaction volume seed
                    this->boundarySeeds_.emplace_back(Helper::makeBoundaryInteractionVolumeSeed(this->problem_(), element, fvGeometry, scvf));

                    // update the index map entries for the global scv faces in the interaction volume
                    for (auto scvfIdxGlobal : this->boundarySeeds_.back().globalScvfIndices())
                    {
                        this->scvfIndexMap_[scvfIdxGlobal] = boundarySeedIndex;
                        isFaceHandled[scvfIdxGlobal] = true;
                    }

                    // increment counter
                    boundarySeedIndex++;
                }
                else
                {
                    // make the inner interaction volume seed only if we are on highest level of all connected elements
                    if (isLocalMaxLevel(element, scvf))
                    {
                        this->seeds_.emplace_back(Helper::makeInnerInteractionVolumeSeed(this->problem_(), element, fvGeometry, scvf));

                        // update the index map entries for the global scv faces in the interaction volume
                        for (auto scvfIdxGlobal : this->seeds_.back().globalScvfIndices())
                        {
                            this->scvfIndexMap_[scvfIdxGlobal] = seedIndex;
                            isFaceHandled[scvfIdxGlobal] = true;
                        }

                        // increment counter
                        seedIndex++;
                    }
                }
            }
        }

        // shrink vectors to actual size
        this->seeds_.shrink_to_fit();
        this->boundarySeeds_.shrink_to_fit();
    }

    bool isLocalMaxLevel(const Element& element, const SubControlVolumeFace& scvf) const
    {
        auto inLevel = element.level();
        for (auto outsideIdx : scvf.outsideScvIndices())
            if (this->problem_().model().globalFvGeometry().element(outsideIdx).level() > inLevel)
                return false;
        return true;
    }
};

} // end namespace


#endif
