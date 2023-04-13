// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkFlux
 * \brief This file contains the data which is required to calculate
 *        diffusive heat fluxes with Fourier's law.
 */
#ifndef DUMUX_FLUX_PNM_FOURIERS_LAW_HH
#define DUMUX_FLUX_PNM_FOURIERS_LAW_HH

#include <dumux/common/math.hh>
#include <dumux/flux/referencesystemformulation.hh>
#include <type_traits>

namespace Dumux::PoreNetwork {

namespace Detail {

struct NoDiffusionType {};

} // end namespace Detail

 /*!
  * \ingroup PoreNetworkFlux
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
