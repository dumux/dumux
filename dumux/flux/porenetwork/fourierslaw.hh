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
 * \ingroup PoreNetworkModels
 * \brief This file contains the data which is required to calculate
 *        diffusive heat fluxes with Fourier's law.
 */
#ifndef DUMUX_FLUX_PNM_FOURIERS_LAW_HH
#define DUMUX_FLUX_PNM_FOURIERS_LAW_HH

#include <dumux/common/math.hh>
#include <type_traits>

namespace Dumux::PoreNetwork {

namespace Detail {

struct NoDiffusionType {};

} // end namespace Detail

 /*!
  * \ingroup PoreNetworkModels
  * \brief Specialization of Fourier's Law for the pore-network model.
  */
template<class MolecularDiffusionType = Detail::NoDiffusionType>
struct PNMFouriersLaw
{

    template<class Problem, class Element, class FVElementGeometry,
             class ElementVolumeVariables, class ElementFluxVariablesCache>
    static auto flux(const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const typename FVElementGeometry::SubControlVolumeFace& scvf,
                     const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        static constexpr auto numPhases = ElementVolumeVariables::VolumeVariables::numFluidPhases();
        using Scalar = typename ElementVolumeVariables::VolumeVariables::PrimaryVariables::value_type;
        Scalar heatflux = 0;

        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            auto insideThermalConducitivity = insideVolVars.fluidThermalConductivity(phaseIdx);
            auto outsideThermalConducitivity = outsideVolVars.fluidThermalConductivity(phaseIdx);

            auto thermalConductivity = Dumux::harmonicMean(insideThermalConducitivity, outsideThermalConducitivity);
            auto area = fluxVarsCache.throatCrossSectionalArea(phaseIdx);

            // calculate the temperature gradient
            const Scalar deltaT = insideVolVars.temperature() - outsideVolVars.temperature();
            const Scalar gradT = deltaT/fluxVarsCache.throatLength();

            heatflux += thermalConductivity*gradT*area;

            if constexpr (!std::is_same_v<MolecularDiffusionType, Detail::NoDiffusionType>)
                heatflux += componentEnthalpyFlux_(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache, phaseIdx);
        }

        return heatflux;
    }

private:
    template<class Problem, class Element, class FVElementGeometry,
             class ElementVolumeVariables, class ElementFluxVariablesCache>
    static auto componentEnthalpyFlux_(const Problem& problem,
                                       const Element& element,
                                       const FVElementGeometry& fvGeometry,
                                       const ElementVolumeVariables& elemVolVars,
                                       const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                       const ElementFluxVariablesCache& elemFluxVarsCache,
                                       const int phaseIdx)
    {
        using Scalar = typename ElementVolumeVariables::VolumeVariables::PrimaryVariables::value_type;
        Scalar heatflux = 0.0;
        using FluidSystem = typename ElementVolumeVariables::VolumeVariables::FluidSystem;
        const auto diffusiveFlux = MolecularDiffusionType::flux(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, elemFluxVarsCache);
        for (int compIdx = 0; compIdx < ElementVolumeVariables::VolumeVariables::numFluidComponents(); ++compIdx)
        {
            const bool insideIsUpstream = diffusiveFlux[compIdx] > 0.0;
            const auto& upstreamVolVars = insideIsUpstream ? elemVolVars[scvf.insideScvIdx()] : elemVolVars[scvf.outsideScvIdx()];
            const Scalar componentEnthalpy = FluidSystem::componentEnthalpy(upstreamVolVars.fluidState(), phaseIdx, compIdx);

            if (MolecularDiffusionType::referenceSystemFormulation() == ReferenceSystemFormulation::massAveraged)
                heatflux += diffusiveFlux[compIdx] * componentEnthalpy;
            else
                heatflux += diffusiveFlux[compIdx] * FluidSystem::molarMass(compIdx) * componentEnthalpy;
        }

        return heatflux;
    }
};

} // end namespace Dumux::PoreNetwork

#endif
