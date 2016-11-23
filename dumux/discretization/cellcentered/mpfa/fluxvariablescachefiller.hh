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
 * \brief The global object of flux var caches
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_GLOBAL_FLUXVARSCACHE_FILLER_HH
#define DUMUX_DISCRETIZATION_CCMPFA_GLOBAL_FLUXVARSCACHE_FILLER_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{
//! Forward declaration of the actual implementation
template<class TypeTag, bool advection, bool diffusion, bool energy>
class CCMpfaFluxVariablesCacheFillerImplementation;

/*!
 * \ingroup ImplicitModel
 * \brief Helper class to fill the flux var caches
 */
template<class TypeTag>
using CCMpfaFluxVariablesCacheFiller = CCMpfaFluxVariablesCacheFillerImplementation<TypeTag,
                                                                                    GET_PROP_VALUE(TypeTag, EnableAdvection),
                                                                                    GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion),
                                                                                    GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;

//! Implementation for only advection
template<class TypeTag>
class CCMpfaFluxVariablesCacheFillerImplementation<TypeTag, true, false, false>
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

    static const int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static const bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static const bool solDependentParams = GET_PROP_VALUE(TypeTag, SolutionDependentParameters);

public:
    //! function to fill the flux var caches
    template<class FluxVarsCacheVector>
    static void fillFluxVarCache(const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const SubControlVolumeFace& scvf,
                                 FluxVarsCacheVector& fluxVarsCache,
                                 const bool updateNeumannOnly = false)
    {
        const bool boundary = problem.model().globalFvGeometry().scvfTouchesBoundary(scvf);

        // if we only want to update the neumann fluxes, skip the rest if we don't touch the boundary
        if (updateNeumannOnly && !boundary)
            return;

        // lambda function to get the permeability tensor
        const auto* prob = &problem;
        auto permFunction = [prob](const Element& element,
                                   const VolumeVariables& volVars,
                                   const SubControlVolume& scv)
                            { return prob->spatialParams().intrinsicPermeability(scv, volVars); };

        // update flux var caches
        if (boundary)
        {
            const auto& seed = problem.model().globalFvGeometry().boundaryInteractionVolumeSeed(scvf);
            BoundaryInteractionVolume iv(seed, problem, fvGeometry, elemVolVars);
            iv.solveLocalSystem(permFunction);

            // we assume phaseIdx = eqIdx here for purely advective problems
            for (unsigned int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            {
                // lambda function defining the upwind factor of the advective flux
                auto advectionUpwindFunction = [eqIdx](const VolumeVariables& volVars)
                                               { return volVars.density(eqIdx)/volVars.viscosity(eqIdx); };
                iv.assembleNeumannFluxes(advectionUpwindFunction, eqIdx);

                // update flux variables cache for the scvf
                const auto scvfIdx = scvf.index();
                auto& cache = fluxVarsCache[scvfIdx];
                cache.updatePhaseNeumannFlux(problem, element, fvGeometry, elemVolVars, scvf, iv, eqIdx);
                if (eqIdx == numEq-1)
                {
                    if (!updateNeumannOnly)
                        cache.updateAdvection(problem, element, fvGeometry, elemVolVars, scvf, iv);
                    cache.setUpdateStatus(true);
                }

                // update flux variable caches of the other scvfs of the interaction volume
                for (const auto scvfIdxJ : iv.globalScvfs())
                {
                    if (scvfIdxJ != scvfIdx)
                    {
                        const auto& scvfJ = fvGeometry.scvf(scvfIdxJ);
                        const auto elementJ = problem.model().globalFvGeometry().element(scvfJ.insideScvIdx());
                        auto& cacheJ = fluxVarsCache[scvfIdxJ];
                        cacheJ.updatePhaseNeumannFlux(problem, elementJ, fvGeometry, elemVolVars, scvfJ, iv, eqIdx);
                        if (eqIdx == numEq-1)
                        {
                            if (!updateNeumannOnly)
                                cacheJ.updateAdvection(problem, elementJ, fvGeometry, elemVolVars, scvfJ, iv);
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
            iv.solveLocalSystem(permFunction);

            // update flux variables cache
            const auto scvfIdx = scvf.index();
            auto& cache = fluxVarsCache[scvfIdx];
            cache.updateAdvection(problem, element, fvGeometry, elemVolVars, scvf, iv);
            cache.setUpdateStatus(true);

            // update flux variable caches of the other scvfs of the interaction volume
            for (const auto scvfIdxJ : iv.globalScvfs())
            {
                if (scvfIdxJ != scvfIdx)
                {
                    const auto& scvfJ = fvGeometry.scvf(scvfIdxJ);
                    const auto elementJ = problem.model().globalFvGeometry().element(scvfJ.insideScvIdx());
                    auto& cacheJ = fluxVarsCache[scvfIdxJ];
                    cacheJ.updateAdvection(problem, elementJ, fvGeometry, elemVolVars, scvfJ, iv);
                    cacheJ.setUpdateStatus(true);
                }
            }
        }
    }

    //! function to update the flux var caches during derivative calculation
    template<class FluxVarsCacheVector>
    static void updateFluxVarCache(const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   const ElementVolumeVariables& elemVolVars,
                                   const SubControlVolumeFace& scvf,
                                   FluxVarsCacheVector& fluxVarsCache)
    {
        if (!solDependentParams)
        {
            fluxVarsCache[scvf.index()].setUpdateStatus(true);

            // if we do not use tpfa on the boundary, we have to update the neumann fluxes
            if (!useTpfaBoundary)
                fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCache, true);
        }
        else
            fillFluxVarCache(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCache);
    }
};

} // end namespace

#endif
