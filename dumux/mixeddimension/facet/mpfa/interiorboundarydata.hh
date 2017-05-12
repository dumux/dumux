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
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FACET_INTERIORBOUNDARYDATA_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FACET_INTERIORBOUNDARYDATA_HH

#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>
#include <dumux/mixeddimension/properties.hh>

namespace Dumux
{

template<class TypeTag>
class CCMpfaFacetCouplingInteriorBoundaryData
{
    // types associated with the low dim domain
    using GlobalProblemTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(GlobalProblemTypeTag, LowDimProblemTypeTag);

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);

    using IndexType = typename GridView::IndexSet::IndexType;
    using LocalIndexType = typename BoundaryInteractionVolume::Traits::LocalIndexType;

    using LowDimProblem = typename GET_PROP_TYPE(LowDimProblemTypeTag, Problem);
    using LowDimSpatialParams = typename GET_PROP_TYPE(LowDimProblemTypeTag, SpatialParams);
    using LowDimVolumeVariables = typename GET_PROP_TYPE(LowDimProblemTypeTag, VolumeVariables);
    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);
    using LowDimElement = typename LowDimGridView::template Codim<0>::Entity;
    using LowDimFVElementGeometry = typename GET_PROP_TYPE(LowDimProblemTypeTag, FVElementGeometry);
    using LowDimSubControlVolume = typename GET_PROP_TYPE(LowDimProblemTypeTag, SubControlVolume);

    //! Dummy type for the CompleteCoupledFacetData struct.
    //! Implementations need to have at least the provided interfaces.
    //! Note that the return types are also "wrong" (Here just to satisfy the compiler)
    struct CompleteCoupledFacetData
    {
        const LowDimProblem& problem() const
        { return *lowDimProblemPtr_; }

        const LowDimSpatialParams& spatialParams() const
        { return *lowDimSpatialParamsPtr_; }

        const LowDimVolumeVariables& volVars() const
        { return lowDimVolVars_; }

        const LowDimElement& element() const
        { return lowDimElement_;}

        const LowDimFVElementGeometry& fvGeometry() const
        { return lowDimFvGeometry_; }

        const LowDimSubControlVolume& scv() const
        { return lowDimScv_; }

        CompleteCoupledFacetData(const LowDimProblem& problem,
                                 const IndexType elementIndex,
                                 LowDimElement&& element,
                                 LowDimFVElementGeometry&& fvGeometry,
                                 LowDimVolumeVariables&& volVars)
        : lowDimProblemPtr_(&problem),
          lowDimSpatialParamsPtr_(&problem.spatialParams()),
          lowDimElement_(std::move(element)),
          lowDimFvGeometry_(std::move(fvGeometry)),
          lowDimVolVars_(std::move(volVars)),
          lowDimScv_(lowDimFvGeometry_.scv(elementIndex))
        {}

    private:
        const LowDimProblem* lowDimProblemPtr_;
        const LowDimSpatialParams* lowDimSpatialParamsPtr_;

        LowDimElement lowDimElement_;
        LowDimFVElementGeometry lowDimFvGeometry_;
        LowDimVolumeVariables lowDimVolVars_;
        const LowDimSubControlVolume& lowDimScv_;
    };

public:
    //! the constructor
    CCMpfaFacetCouplingInteriorBoundaryData(const Problem& problem,
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
    LowDimVolumeVariables facetVolVars(const FVElementGeometry& fvGeometry) const
    {
        return problem_().couplingManager().lowDimVolVars(problem_().model().globalFvGeometry().element(elementIndex()),
                                                          fvGeometry,
                                                          fvGeometry.scvf(scvfIndex()));
    }

    //! returns the volume variables for interior dirichlet boundaries
    LowDimVolumeVariables facetVolVars(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf) const
    {
        assert(scvf.index() == scvfIndex() && "calling facet volume variables for an scvf other than the bound one");
        return problem_().couplingManager().lowDimVolVars(problem_().model().globalFvGeometry().element(elementIndex()),
                                                          fvGeometry,
                                                          scvf);
    }

    //! The following interface is here for compatibility reasons to be overloaded for problems using facet coupling.
    //! prepares all the necessary variables of the other domain.
    //! Note that also an implementation of the CompleteFacetData structure has to be provided.
    CompleteCoupledFacetData completeCoupledFacetData(const FVElementGeometry& fvGeometry) const
    {
        const auto& couplingMapper = problem_().couplingManager().couplingMapper();

        // get coupling data for this scvf
        const auto element = problem_().model().globalFvGeometry().element(elementIndex());
        const auto lowDimElementIndex = couplingMapper.getBulkCouplingData(element).getCouplingLowDimElementIndex(scvfIndex());

        // obtain data necessary to fully instantiate the complete coupled facet data
        const auto& lowDimProblem = problem_().couplingManager().lowDimProblem();
        auto lowDimElement = lowDimProblem.model().globalFvGeometry().element(lowDimElementIndex);

        auto lowDimFvGeometry = localView(lowDimProblem.model().globalFvGeometry());
        lowDimFvGeometry.bindElement(lowDimElement);

        LowDimVolumeVariables lowDimVolVars;
        lowDimVolVars.update(lowDimProblem.model().elementSolution(lowDimElement, lowDimProblem.model().curSol()),
                             lowDimProblem,
                             lowDimElement,
                             lowDimFvGeometry.scv(lowDimElementIndex));

        return CompleteCoupledFacetData(lowDimProblem,
                                        lowDimElementIndex,
                                        std::move(lowDimElement),
                                        std::move(lowDimFvGeometry),
                                        std::move(lowDimVolVars));
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
