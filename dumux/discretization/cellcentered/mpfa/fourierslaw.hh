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
 *        heat conduction fluxes with Fourier's law for cell-centered MPFA models.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FOURIERS_LAW_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

/*!
 * \ingroup FouriersLaw
 * \brief Specialization of Fourier's Law for the CCMpfa method.
 */
template <class TypeTag>
class FouriersLawImplementation<TypeTag, DiscretizationMethods::CCMpfa>
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVarsCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);

    // Always use the dynamic type for vectors (compatibility with the boundary)
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using CoefficientVector = typename BoundaryInteractionVolume::Vector;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;

    static constexpr bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static constexpr bool enableInteriorBoundaries = GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries);

    static constexpr int energyEqIdx = GET_PROP_TYPE(TypeTag, Indices)::energyEqIdx;

    //! The cache used in conjunction with the mpfa Fourier's Law
    class MpfaFouriersLawCache
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

        //! If the useTpfaBoundary property is set to false, the boundary conditions
        //! are put into the local systems leading to possible contributions on all faces
        Scalar heatNeumannFlux() const
        { return heatNeumannFlux_; }

    private:
        // Quantities associated with heat conduction
        Stencil heatConductionVolVarsStencil_;
        CoefficientVector heatConductionTij_;
        Scalar heatNeumannFlux_;
    };

    //! Class that fills the cache corresponding to mpfa Darcy's Law
    class MpfaFouriersLawCacheFiller
    {
    public:
        //! Function to fill an MpfaDarcysLawCache of a given scvf
        //! This interface has to be met by any cache filler class for heat conduction quantities
        template<class FluxVariablesCacheFiller>
        static void fill(FluxVariablesCache& scvfFluxVarsCache,
                         const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace& scvf,
                         const FluxVariablesCacheFiller& fluxVarsCacheFiller)
        {
            // get interaction volume from the flux vars cache filler & upate the cache
            if (problem.model().globalFvGeometry().isInBoundaryInteractionVolume(scvf))
                scvfFluxVarsCache.updateHeatConduction(fluxVarsCacheFiller.boundaryInteractionVolume(), scvf);
            else
                scvfFluxVarsCache.updateHeatConduction(fluxVarsCacheFiller.interactionVolume(), scvf);
        }
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCMpfa;

    // state the type for the corresponding cache and its filler
    using Cache = MpfaFouriersLawCache;
    using CacheFiller = MpfaFouriersLawCacheFiller;

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

        const bool isInteriorBoundary = enableInteriorBoundaries && fluxVarsCache.isInteriorBoundary();
        // For interior Neumann boundaries when using Tpfa on boundaries, return the user-specified flux
        if (isInteriorBoundary
            && useTpfaBoundary
            && fluxVarsCache.interiorBoundaryDataSelf().faceType() == MpfaFaceTypes::interiorNeumann)
            return scvf.area()*
                   elemVolVars[scvf.insideScvIdx()].extrusionFactor()*
                   problem.neumann(element, fvGeometry, elemVolVars, scvf)[energyEqIdx];

        // calculate Tij*tj
        Scalar flux(0.0);
        unsigned int localIdx = 0;
        for (const auto volVarIdx : volVarsStencil)
            flux += tij[localIdx++]*elemVolVars[volVarIdx].temperature();

        // if no interior boundaries are present, return heat conduction flux
        if (!enableInteriorBoundaries)
            return useTpfaBoundary ? flux : flux + fluxVarsCache.heatNeumannFlux();

        // Handle interior boundaries
        flux += Implementation::computeInteriorBoundaryContribution(element, fvGeometry, fluxVarsCache);

        // return overall resulting flux
        return useTpfaBoundary ? flux : flux + fluxVarsCache.heatNeumannFlux();
    }

    static Scalar computeInteriorBoundaryContribution(const Element& element,
                                                      const FVElementGeometry& fvGeometry,
                                                      const FluxVariablesCache& fluxVarsCache)
    {
        // obtain the transmissibilites associated with all pressures
        const auto& tij = fluxVarsCache.heatConductionTij();

        // the interior dirichlet boundaries local indices start after
        // the cell and the domain Dirichlet boundary pressures
        const auto startIdx = fluxVarsCache.heatConductionVolVarsStencil().size();

        // add interior Dirichlet boundary contributions
        Scalar flux = 0.0;
        for (auto&& data : fluxVarsCache.interiorBoundaryData())
            if (data.faceType() == MpfaFaceTypes::interiorDirichlet)
                flux += tij[startIdx + data.localIndexInInteractionVolume()]*data.facetVolVars(element, fvGeometry).temperature();

        return flux;
    }
};

} // end namespace Dumux

#endif
