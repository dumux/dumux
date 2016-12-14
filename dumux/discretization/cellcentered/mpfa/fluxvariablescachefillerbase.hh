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
 * \brief Base classes for the flux variables cache filler
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_GLOBAL_FLUXVARSCACHE_FILLER_BASE_HH
#define DUMUX_DISCRETIZATION_CCMPFA_GLOBAL_FLUXVARSCACHE_FILLER_BASE_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

//! Forward declaration of property tags
namespace Properties
{
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(NumComponents);
NEW_PROP_TAG(EffectiveDiffusivityModel);
}

//! Fills the advective quantities in the caches
template<class TypeTag>
class CCMpfaAdvectionCacheFiller
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using Element = typename GridView::template Codim<0>::Entity;

    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static const bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static const bool enableDiffusion = GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion);

public:
    //! function to fill the flux var caches
    template<class FluxVarsCacheContainer>
    static void fillCaches(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf,
                           FluxVarsCacheContainer& fluxVarsCacheContainer)
    {
        // lambda function to get the permeability tensor
        auto getK = [&problem](const Element& element,
                               const VolumeVariables& volVars,
                               const SubControlVolume& scv)
        { return problem.spatialParams().intrinsicPermeability(scv, volVars); };

        // update flux var caches
        if (problem.model().globalFvGeometry().scvfTouchesBoundary(scvf))
        {
            const auto& seed = problem.model().globalFvGeometry().boundaryInteractionVolumeSeed(scvf);
            BoundaryInteractionVolume iv(seed, problem, fvGeometry, elemVolVars);
            iv.solveLocalSystem(getK);

            // update the transmissibilities etc for each phase
            const auto scvfIdx = scvf.index();
            auto& cache = fluxVarsCacheContainer[scvfIdx];
            cache.updateAdvection(scvf, iv);
            cache.setUpdateStatus(true);

            //! set the neumann boundary conditions in case we do not use tpfa on the
            //! boundary and diffusion is not enabled (then we assume neumann BCs to be diffusive)
            if (!useTpfaBoundary && !enableDiffusion)
                for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                {
                    iv.assembleNeumannFluxes(phaseIdx);
                    cache.updatePhaseNeumannFlux(scvf, iv, phaseIdx);
                }

            for (const auto scvfIdxJ : iv.globalScvfs())
            {
                if (scvfIdxJ != scvfIdx)
                {
                    const auto& scvfJ = fvGeometry.scvf(scvfIdxJ);
                    auto& cacheJ = fluxVarsCacheContainer[scvfIdxJ];
                    cacheJ.updateAdvection(scvfJ, iv);
                    cacheJ.setUpdateStatus(true);

                    if (!useTpfaBoundary && !enableDiffusion)
                        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                        {
                            iv.assembleNeumannFluxes(phaseIdx);
                            cacheJ.updatePhaseNeumannFlux(scvfJ, iv, phaseIdx);
                        }
                }
            }
        }
        else
        {
            const auto& seed = problem.model().globalFvGeometry().interactionVolumeSeed(scvf);
            InteractionVolume iv(seed, problem, fvGeometry, elemVolVars);
            iv.solveLocalSystem(getK);

            // update flux variables cache
            const auto scvfIdx = scvf.index();
            auto& cache = fluxVarsCacheContainer[scvfIdx];
            cache.updateAdvection(scvf, iv);
            cache.setUpdateStatus(true);

            // update flux variable caches of the other scvfs of the interaction volume
            for (const auto scvfIdxJ : iv.globalScvfs())
            {
                if (scvfIdxJ != scvfIdx)
                {
                    const auto& scvfJ = fvGeometry.scvf(scvfIdxJ);
                    auto& cacheJ = fluxVarsCacheContainer[scvfIdxJ];
                    cacheJ.updateAdvection(scvfJ, iv);
                    cacheJ.setUpdateStatus(true);
                }
            }
        }
    }
};

//! Fills the diffusive quantities in the caches
template<class TypeTag>
class CCMpfaDiffusionCacheFiller
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

    using Element = typename GridView::template Codim<0>::Entity;

    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static const int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);

public:
    //! function to fill the flux var caches
    template<class FluxVarsCacheContainer>
    static void fillCaches(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf,
                           FluxVarsCacheContainer& fluxVarsCacheContainer)
    {
        if (problem.model().globalFvGeometry().scvfTouchesBoundary(scvf))
        {
            const auto& seed = problem.model().globalFvGeometry().boundaryInteractionVolumeSeed(scvf);
            BoundaryInteractionVolume iv(seed, problem, fvGeometry, elemVolVars);

            //! update the cache for all components in all phases. We exclude the case
            //! phaseIdx = compIdx here, as diffusive fluxes of the major component in its
            //! own phase are not calculated explicitly during assembly (see compositional local residual)
            for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                for (unsigned int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    if (phaseIdx == compIdx)
                        continue;

                    // lambda function to get the diffusion coefficient/tensor
                    auto getD = [phaseIdx, compIdx](const Element& element,
                                                    const VolumeVariables& volVars,
                                                    const SubControlVolume& scv)
                    { return volVars.diffusionCoefficient(phaseIdx, compIdx); };

                    // solve for the transmissibilities
                    iv.solveLocalSystem(getD);

                    // update the caches
                    const auto scvfIdx = scvf.index();
                    auto& cache = fluxVarsCacheContainer[scvfIdx];
                    cache.updateDiffusion(scvf, iv, phaseIdx, compIdx);
                    cache.setUpdateStatus(true);

                    for (const auto scvfIdxJ : iv.globalScvfs())
                    {
                        if (scvfIdxJ != scvfIdx)
                        {
                            const auto& scvfJ = fvGeometry.scvf(scvfIdxJ);
                            auto& cacheJ = fluxVarsCacheContainer[scvfIdxJ];
                            cacheJ.updateDiffusion(scvfJ, iv, phaseIdx, compIdx);
                            cacheJ.setUpdateStatus(true);
                        }
                    }
                }
            }
        }
        else
        {
            const auto& seed = problem.model().globalFvGeometry().interactionVolumeSeed(scvf);
            InteractionVolume iv(seed, problem, fvGeometry, elemVolVars);

            //! update the cache for all components in all phases. We exclude the case
            //! phaseIdx = compIdx here, as diffusive fluxes of the major component in its
            //! own phase are not calculated explicitly during assembly (see compositional local residual)
            for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                for (unsigned int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    if (phaseIdx == compIdx)
                        continue;

                    // lambda function to get the diffusion coefficient/tensor
                    // TODO: How to include effective diffusivity?
                    auto getD = [phaseIdx, compIdx](const Element& element,
                                                    const VolumeVariables& volVars,
                                                    const SubControlVolume& scv)
                    { return volVars.diffusionCoefficient(phaseIdx, compIdx); };

                    // solve for the transmissibilities
                    iv.solveLocalSystem(getD);

                    // update the caches
                    const auto scvfIdx = scvf.index();
                    auto& cache = fluxVarsCacheContainer[scvfIdx];
                    cache.updateDiffusion(scvf, iv, phaseIdx, compIdx);
                    cache.setUpdateStatus(true);

                    for (const auto scvfIdxJ : iv.globalScvfs())
                    {
                        if (scvfIdxJ != scvfIdx)
                        {
                            const auto& scvfJ = fvGeometry.scvf(scvfIdxJ);
                            auto& cacheJ = fluxVarsCacheContainer[scvfIdxJ];
                            cacheJ.updateDiffusion(scvfJ, iv, phaseIdx, compIdx);
                            cacheJ.setUpdateStatus(true);
                        }
                    }
                }
            }
        }
    }
};



} // end namespace

#endif
