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
 * \brief A class to store info on interior boundaries
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_INTERIORBOUNDARYDATA_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_INTERIORBOUNDARYDATA_HH

#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>

namespace Dumux
{
/*!
 * \ingroup CCMpfa
 * \brief Class to store additional data related to interior boundary facets.
 */
template<class TypeTag>
class InteriorBoundaryData
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);

    using IndexType = typename GridView::IndexSet::IndexType;
    using LocalIndexType = typename BoundaryInteractionVolume::LocalIndexType;

public:
    //! the constructor
    InteriorBoundaryData(const Problem& problem,
                         IndexType elementIndex,
                         IndexType scvfIndex,
                         LocalIndexType localIndex,
                         MpfaFaceTypes faceType)
    : problemPtr_(&problem),
      elementIndex_(elementIndex),
      scvfIndex_(scvfIndex),
      localIndex_(localIndex),
      faceType_(faceType)
    {}

    //! returns the global index of the element/scv connected to the interior boundary
    IndexType elementIndex() const
    { return elementIndex_; }

    //! returns the global index of the scvf connected to the interior boundary
    IndexType scvfIndex() const
    { return scvfIndex_; }

    //! returns the local index i of the scvf within the interaction volume.
    //! This is either:
    //!   - the i-th flux face index (interior neumann boundaries)
    //!   - the i-th interior dirichlet face (interior dirichlet boundaries)
    LocalIndexType localIndexInInteractionVolume() const
    { return localIndex_; }

    //! returns the face type of this scvf
    MpfaFaceTypes faceType() const
    { return faceType_; }

    //! returns the volume variables for interior dirichlet boundaries
    VolumeVariables facetVolVars(const Element& element, const FVElementGeometry& fvGeometry) const
    {
        //! This can only be called for interior Dirichlet boundaries
        assert(faceType_ == MpfaFaceTypes::interiorDirichlet && "requesting Dirichlet vol vars for a face which is"
                                                                "not marked as interior Dirichlet face.");

        auto priVars = problem_().dirichlet(element, fvGeometry.scvf(scvfIndex()));

        VolumeVariables volVars;
        volVars.update(ElementSolutionVector({priVars}),
                       problem_(),
                       element,
                       fvGeometry.scv(elementIndex()));

        return volVars;
    }

private:
    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;

    IndexType elementIndex_;
    IndexType scvfIndex_;
    LocalIndexType localIndex_;
    MpfaFaceTypes faceType_;
};
} // end namespace

#endif
