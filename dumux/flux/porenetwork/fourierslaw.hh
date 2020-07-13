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
 *        diffusive heat fluxes with Fourier's law.
 */
#ifndef DUMUX_FLUX_PNM_FOURIERS_LAW_HH
#define DUMUX_FLUX_PNM_FOURIERS_LAW_HH

#include <dumux/common/math.hh>

namespace Dumux::PoreNetwork {

 /*!
  * \ingroup PoreNetworkFlux
  * \brief Specialization of Fourier's Law for the pore-network model.
  */
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
        }
        return heatflux;
    }
};

} // end namespace Dumux::PoreNetwork

#endif
