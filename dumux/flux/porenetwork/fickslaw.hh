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

            auto insideDiffCoeff = getDiffusionCoefficient_(phaseIdx, compIdx, insideVolVars);
            auto outsideDiffCoeff = getDiffusionCoefficient_(phaseIdx, compIdx, outsideVolVars);

            // scale by extrusion factor
            insideDiffCoeff *= insideVolVars.extrusionFactor();
            outsideDiffCoeff *= outsideVolVars.extrusionFactor();

            // the resulting averaged diffusion coefficient
            const auto diffCoeff = harmonicMean(insideDiffCoeff, outsideDiffCoeff);

            const Scalar insideMoleFraction = massOrMoleFraction(insideVolVars, referenceSystem, phaseIdx, compIdx);
            const Scalar outsideMoleFraction = massOrMoleFraction(outsideVolVars, referenceSystem, phaseIdx, compIdx);

            componentFlux[compIdx] = density * (insideMoleFraction - outsideMoleFraction) / throatLength * diffCoeff * phaseCrossSectionalArea;
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
