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
 *        diffusive mass fluxes due to molecular diffusion with Fick's law.
 */
#ifndef DUMUX_FLUX_PNM_FICKS_LAW_HH
#define DUMUX_FLUX_PNM_FICKS_LAW_HH

#include <dune/common/fvector.hh>
#include <dumux/common/math.hh>
#include <dumux/flux/referencesystemformulation.hh>
#include <dumux/flux/fickiandiffusioncoefficients.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PoreNetworkFlux
 * \brief Specialization of Fick's Law for the pore-network model.
 */
template<class Scalar, int numPhases, int numComponents,
         ReferenceSystemFormulation referenceSystem = ReferenceSystemFormulation::massAveraged>
class PNMFicksLaw
{
public:
    using DiffusionCoefficientsContainer = FickianDiffusionCoefficients<Scalar, numPhases, numComponents>;

    //return the reference system
    static constexpr ReferenceSystemFormulation referenceSystemFormulation()
    { return referenceSystem; }

    template<class Problem, class Element, class FVElementGeometry,
             class ElementVolumeVariables, class ElementFluxVariablesCache>
    static Dune::FieldVector<Scalar, numComponents>
    flux(const Problem& problem,
         const Element& element,
         const FVElementGeometry& fvGeometry,
         const ElementVolumeVariables& elemVolVars,
         const typename FVElementGeometry::SubControlVolumeFace& scvf,
         const int phaseIdx,
         const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        Dune::FieldVector<Scalar, numComponents> componentFlux(0.0);

        // get inside and outside diffusion tensors and calculate the harmonic mean
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        const Scalar density = 0.5 * (massOrMolarDensity(insideVolVars, referenceSystem, phaseIdx) + massOrMolarDensity(outsideVolVars, referenceSystem, phaseIdx));
        const Scalar throatLength = fluxVarsCache.throatLength();
        const Scalar phaseCrossSectionalArea = fluxVarsCache.throatCrossSectionalArea(phaseIdx);

        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            if(compIdx == phaseIdx)
                continue;

            auto insideD = getDiffusionCoefficient_(phaseIdx, compIdx, insideVolVars);
            auto outsideD = getDiffusionCoefficient_(phaseIdx, compIdx, outsideVolVars);

            // scale by extrusion factor
            insideD *= insideVolVars.extrusionFactor();
            outsideD *= outsideVolVars.extrusionFactor();

            // the resulting averaged diffusion coefficient
            const auto D = harmonicMean(insideD, outsideD);

            const Scalar insideMoleFraction = massOrMoleFraction(insideVolVars, referenceSystem, phaseIdx, compIdx);
            const Scalar outsideMoleFraction = massOrMoleFraction(outsideVolVars, referenceSystem, phaseIdx, compIdx);

            componentFlux[compIdx] = density * (insideMoleFraction - outsideMoleFraction) / throatLength * D * phaseCrossSectionalArea;
            componentFlux[phaseIdx] -= componentFlux[compIdx];
        }
        return componentFlux;
    }
private:

    template<class VolumeVariables>
    static Scalar getDiffusionCoefficient_(const int phaseIdx, const int compIdx,
                                           const VolumeVariables& volVars)
    {
        using FluidSystem = typename VolumeVariables::FluidSystem;

        if constexpr (!FluidSystem::isTracerFluidSystem())
        {
            const auto mainCompIdx = FluidSystem::getMainComponent(phaseIdx);
            return volVars.diffusionCoefficient(phaseIdx, mainCompIdx, compIdx);
        }
        else
            return volVars.diffusionCoefficient(0, 0, compIdx);
    }
};
} // end namespace Dumux::PoreNetwork

#endif
