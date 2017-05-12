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
 * \brief Base classes for interaction volumes of mpfa models with active coupling over the element facets.
 */
#ifndef DUMUX_MIXEDDIMENSION_FACET_MPFA_O_INTERACTIONVOLUME_HH
#define DUMUX_MIXEDDIMENSION_FACET_MPFA_O_INTERACTIONVOLUME_HH

#include <dumux/discretization/cellcentered/mpfa/interactionvolume.hh>
#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

namespace Dumux
{
//! Forward declaration
template<class TypeTag> class CCMpfaOFacetCouplingInteractionVolume;

//! Specialization of the interaction volume traits class for the o-method in coupled models
template<class TypeTag>
class CCMpfaOFacetCouplingInteractionVolumeTraits : public CCMpfaOInteractionVolumeTraits<TypeTag>
{
public:
    using BoundaryInteractionVolume = CCMpfaOFacetCouplingInteractionVolume<TypeTag>;
};

// the o-method interaction volume is substituted by the one including data on the facet element's
// tensorial quantities into the local system to be solved.
template<class TypeTag>
class CCMpfaOFacetCouplingInteractionVolume : public CCMpfaOInteractionVolume<TypeTag,
                                                                              CCMpfaOFacetCouplingInteractionVolumeTraits<TypeTag>,
                                                                              CCMpfaOFacetCouplingInteractionVolume<TypeTag>>
{
    using Traits = CCMpfaOFacetCouplingInteractionVolumeTraits<TypeTag>;
    using ThisType = CCMpfaOFacetCouplingInteractionVolume<TypeTag>;
    using ParentType = CCMpfaOInteractionVolume<TypeTag, Traits, ThisType>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteriorBoundaryData = typename GET_PROP_TYPE(TypeTag, InteriorBoundaryData);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using LocalScvfType = typename Traits::LocalScvfType;
    using LocalIndexType = typename Traits::LocalIndexType;

    static constexpr bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);

public:

    // We include data on the tensorial quantities of the facet elements here
    template<typename GetTensorFunction>
    Scalar interiorNeumannTerm(const GetTensorFunction& getTensor,
                               const Element& element,
                               const LocalScvfType& localScvf,
                               const InteriorBoundaryData& data) const
    {
        // obtain the complete data on the facet element
        const auto completeFacetData = data.completeCoupledFacetData(this->fvGeometry());

        // calculate "leakage factor"
        const auto n = localScvf.unitOuterNormal();
        const auto v = [&] ()
                        {
                            auto res = n;
                            res *= -0.5*completeFacetData.volVars().extrusionFactor();
                            res /= res.two_norm2();
                            return res;
                        } ();

        // substract (n*T*v)*Area from diagonal matrix entry
        const auto facetTensor = getTensor(completeFacetData.problem(),
                                           completeFacetData.element(),
                                           completeFacetData.volVars(),
                                           completeFacetData.fvGeometry(),
                                           completeFacetData.scv());

        return localScvf.area()*
               this->elemVolVars()[localScvf.insideGlobalScvIndex()].extrusionFactor()*
               MpfaHelper::nT_M_v(n, facetTensor, v);
    }

    void assembleNeumannFluxVector()
    {
        // initialize the neumann fluxes vector to zero
        this->neumannFluxes_.resize(this->fluxScvfIndexSet().size(), PrimaryVariables(0.0));

        if (!this->onDomainOrInteriorBoundary() || useTpfaBoundary)
            return;

        LocalIndexType fluxFaceIdx = 0;
        for (auto localFluxFaceIdx : this->fluxScvfIndexSet())
        {
            const auto& localScvf = this->localScvf(localFluxFaceIdx);
            const auto faceType = localScvf.faceType();

            if (faceType == MpfaFaceTypes::neumann)
            {
                const auto& element = this->localElement(localScvf.insideLocalScvIndex());
                const auto& globalScvf = this->fvGeometry().scvf(localScvf.insideGlobalScvfIndex());
                auto neumannFlux = this->problem().neumann(element, this->fvGeometry(), this->elemVolVars(), globalScvf);
                neumannFlux *= globalScvf.area();
                neumannFlux *= this->elemVolVars()[globalScvf.insideScvIdx()].extrusionFactor();

                // The flux is assumed to be prescribed in the form of -D*gradU
                this->neumannFluxes_[fluxFaceIdx] = neumannFlux;
            }

            fluxFaceIdx++;
        }
    }
};

} // end namespace

#endif
