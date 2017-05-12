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
 * \brief This file contains the data which is required to calculate
 *        heat conduction fluxes with Fourier's law for cell-centered MPFA models
 *        in the presence of lower dimensional (coupled) elements living on this
 *        domain's element facets.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FACET_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FACET_FOURIERS_LAW_HH

#include <dumux/discretization/cellcentered/mpfa/fourierslaw.hh>
#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>

namespace Dumux
{

/*!
 * \ingroup FouriersLaw
 * \brief Specialization of Fourier's Law for the CCMpfa method with lower dimensional
 *        elements living on the bulk elements' facets.
 */
template <class TypeTag>
class CCMpfaFacetCouplingFouriersLaw : public FouriersLawImplementation<TypeTag, DiscretizationMethods::CCMpfa>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVarsCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);

    // Always use the dynamic type for vectors (compatibility with the boundary)
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;

    static constexpr bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static constexpr int energyEqIdx = GET_PROP_TYPE(TypeTag, Indices)::energyEqIdx;

    //! The cache used in conjunction with the mpfa Fourier's Law
    class MpfaFacetCouplingFouriersLawCache
    {
        using Stencil = typename BoundaryInteractionVolume::GlobalIndexSet;
    public:
        // update cached objects for heat conduction
        template<typename InteractionVolume>
        void updateHeatConduction(const InteractionVolume& iv, const SubControlVolumeFace &scvf)
        {
            const auto& localFaceData = iv.getLocalFaceData(scvf);
            heatConductionVolVarsStencil_ = iv.volVarsStencil();
            heatConductionTij_ = iv.getTransmissibilities(localFaceData);
            heatNeumannFlux_ = iv.getNeumannFlux(localFaceData, energyEqIdx);
            heatConductionCij_ = iv.getNeumannFluxTransformationCoefficients(localFaceData);
        }

        //! Returns the volume variables indices necessary for heat conduction flux
        //! computation. This includes all participating boundary volume variables
        //! and it can be different for the phases & components.
        const Stencil& heatConductionVolVarsStencil() const
        { return heatConductionVolVarsStencil_; }

        //! Returns the transmissibilities associated with the volume variables
        //! This can be different for the phases & components.
        const CoefficientVector& heatConductionTij() const
        { return heatConductionTij_; }

        //! Returns the vector of coefficients with which the vector of neumann boundary conditions
        //! has to be multiplied in order to transform them on the scvf this cache belongs to
        const CoefficientVector& heatConductionCij() const
        { return heatConductionCij_; }

        //! If the useTpfaBoundary property is set to false, the boundary conditions
        //! are put into the local systems leading to possible contributions on all faces
        Scalar heatNeumannFlux() const
        { return heatNeumannFlux_; }

    private:
        // Quantities associated with heat conduction
        Stencil heatConductionVolVarsStencil_;
        CoefficientVector heatConductionTij_;
        CoefficientVector heatConductionCij_;
        Scalar heatNeumannFlux_;
    };
    using CoefficientVector = typename BoundaryInteractionVolume::Traits::Vector;

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCMpfa;

    // state the new type for the corresponding cache
    using Cache = MpfaFacetCouplingFouriersLawCache;

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto& volVarsStencil = fluxVarsCache.heatConductionVolVarsStencil();
        const auto& tij = fluxVarsCache.heatConductionTij();

        // calculate Tij*tj
        Scalar flux(0.0);
        unsigned int localIdx = 0;
        for (const auto volVarIdx : volVarsStencil)
            flux += tij[localIdx++]*elemVolVars[volVarIdx].temperature();

        // Handle interior boundaries
        flux += computeInteriorBoundaryContribution(fvGeometry, elemVolVars, fluxVarsCache);

        // return overall resulting flux
        return useTpfaBoundary ? flux : flux + fluxVarsCache.heatNeumannFlux();
    }

    static Scalar computeInteriorBoundaryContribution(const FVElementGeometry& fvGeometry,
                                                      const ElementVolumeVariables& elemVolVars,
                                                      const FluxVariablesCache& fluxVarsCache)
    {
        // obtain the transmissibilites associated with all pressures
        const auto& tij = fluxVarsCache.heatConductionTij();

        // the interior dirichlet boundaries local indices start after
        // the cell and the domain Dirichlet boundary pressures
        const auto startIdx = fluxVarsCache.heatConductionVolVarsStencil().size();

        // The vector of interior neumann fluxes
        const auto& cij = fluxVarsCache.heatConductionCij();
        CoefficientVector facetCouplingFluxes(cij.size(), 0.0);

        // add interior Dirichlet boundary contributions
        Scalar flux = 0.0;
        for (auto&& data : fluxVarsCache.interiorBoundaryData())
        {
            // Add additional Dirichlet fluxes for interior Dirichlet faces
            if (data.faceType() == MpfaFaceTypes::interiorDirichlet)
                flux += tij[startIdx + data.localIndexInInteractionVolume()]*data.facetVolVars(fvGeometry).temperature();

            // add neumann contributions
            if (data.faceType() == MpfaFaceTypes::interiorNeumann)
            {
                // get the scvf corresponding to actual interior neumann face
                const auto& curScvf = fvGeometry.scvf(data.scvfIndex());

                // get the complete data of the actual interior neumann face
                const auto completeFacetData = data.completeCoupledFacetData(fvGeometry);

                // calculate "leakage factor"
                const auto n = curScvf.unitOuterNormal();
                const auto v = [&] ()
                                {
                                    auto res = n;
                                    res *= -0.5*completeFacetData.volVars().extrusionFactor();
                                    res /= res.two_norm2();
                                    return res;
                                } ();

                // get the thermal conductivity in the facet element
                const auto facetLambda = ThermalConductivityModel::effectiveThermalConductivity(completeFacetData.volVars(),
                                                                                                completeFacetData.spatialParams(),
                                                                                                completeFacetData.element(),
                                                                                                completeFacetData.fvGeometry(),
                                                                                                completeFacetData.scv());
                // add value to vector of interior neumann fluxes
                facetCouplingFluxes[data.localIndexInInteractionVolume()] += completeFacetData.volVars().temperature()*
                                                                             curScvf.area()*
                                                                             elemVolVars[curScvf.insideScvIdx()].extrusionFactor()*
                                                                             MpfaHelper::nT_M_v(n, facetLambda, v);
            }
        }

        return flux + cij*facetCouplingFluxes;
    }
};

} // end namespace Dumux

#endif
