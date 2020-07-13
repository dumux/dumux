// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 *        the fluxes of the pore network model over a face of a finite volume.
 */
#ifndef DUMUX_FLUX_PNM_ADVECTION_HH
#define DUMUX_FLUX_PNM_ADVECTION_HH

namespace Dumux::PoreNetwork::Detail {

template<class... TransmissibilityLawTypes>
struct Transmissibility : public TransmissibilityLawTypes... {};

} // end namespace Dumux::PoreNetwork::Detail

namespace Dumux::PoreNetwork {

/*!
 * \file
 * \ingroup PoreNetworkFlux
 * \brief Hagen–Poiseuille-type flux law to describe the advective flux for pore-network models.
 */
template<class ScalarT, class... TransmissibilityLawTypes>
class CreepingFlow
{

public:
    //! Export the Scalar type
    using Scalar = ScalarT;

    //! Export the transmissibility law
    using Transmissibility = Detail::Transmissibility<TransmissibilityLawTypes...>;

    /*!
     * \brief Returns the advective flux of a fluid phase
     *        across the given sub-control volume face (corresponding to a pore throat).
     * \note The flux is given in N*m, and can be converted
     *       into a volume flux (m^3/s) or mass flux (kg/s) by applying an upwind scheme
     *       for the mobility (1/viscosity) or the product of density and mobility, respectively.
     */
    template<class Problem, class Element, class FVElementGeometry,
             class ElementVolumeVariables, class SubControlVolumeFace, class ElemFluxVarsCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const int phaseIdx,
                       const ElemFluxVarsCache& elemFluxVarsCache)
    {
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        // Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];

        // calculate the pressure difference
        const Scalar deltaP = insideVolVars.pressure(phaseIdx) - outsideVolVars.pressure(phaseIdx);
        const Scalar transmissibility = fluxVarsCache.transmissibility(phaseIdx);

        Scalar volumeFlow = transmissibility*deltaP;

        // add gravity term
        static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");
        if (enableGravity)
        {
            const Scalar rho = 0.5*insideVolVars.density(phaseIdx) + 0.5*outsideVolVars.density(phaseIdx);
            const Scalar g = problem.spatialParams().gravity(scvf.center()) * scvf.unitOuterNormal();

            // The transmissibility is with respect to the effective throat length (potentially dropping the pore body radii).
            // For gravity, we need to consider the total throat length (i.e., the cell-center to cell-center distance).
            // This might cause some inconsistencies TODO: is there a better way?
            volumeFlow += transmissibility * fluxVarsCache.poreToPoreDistance() * rho * g;
        }

        return volumeFlow;
    }

    /*!
     * \brief Returns the throat conductivity
     */
    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
    static Scalar calculateTransmissibility(const Problem& problem,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                            const ElementVolumeVariables& elemVolVars,
                                            const FluxVariablesCache& fluxVarsCache,
                                            const int phaseIdx)
    {
        if constexpr (ElementVolumeVariables::VolumeVariables::numFluidPhases() == 1)
            return Transmissibility::singlePhaseTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVarsCache, phaseIdx);
        else
        {
            static_assert(ElementVolumeVariables::VolumeVariables::numFluidPhases() == 2);

            const auto& spatialParams = problem.spatialParams();
            using FluidSystem = typename ElementVolumeVariables::VolumeVariables::FluidSystem;
            const int wPhaseIdx = spatialParams.template wettingPhase<FluidSystem>(element, elemVolVars);
            const bool invaded = fluxVarsCache.invaded();

            if (phaseIdx == wPhaseIdx)
            {
                return invaded ? Transmissibility::wettingLayerTransmissibility(element, fvGeometry, scvf, fluxVarsCache)
                               : Transmissibility::singlePhaseTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVarsCache, phaseIdx);
            }
            else // non-wetting phase
            {
                return invaded ? Transmissibility::nonWettingPhaseTransmissibility(element, fvGeometry, scvf, fluxVarsCache)
                               : 0.0;
            }
        }
    }

    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
    static std::array<Scalar, 2> calculateTransmissibilities(const Problem& problem,
                                                             const Element& element,
                                                             const FVElementGeometry& fvGeometry,
                                                             const ElementVolumeVariables& elemVolVars,
                                                             const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                                             const FluxVariablesCache& fluxVarsCache)
    {
        static_assert(ElementVolumeVariables::VolumeVariables::numFluidPhases() == 1);
        const Scalar t = calculateTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVarsCache, 0);
        return {t, -t};
    }
};


} // end namespace Dumux

#endif
