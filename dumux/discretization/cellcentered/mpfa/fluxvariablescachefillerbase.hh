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
#include <dumux/discretization/methods.hh>

namespace Dumux
{

//! Forward declaration of property tags
namespace Properties
{
NEW_PROP_TAG(Indices);
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(NumComponents);
NEW_PROP_TAG(ThermalConductivityModel);
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
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static constexpr bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static constexpr bool enableDiffusion = GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion);

public:
    //! function to fill the flux var caches
    template<class FluxVarsCacheContainer, class InteractionVolume>
    static void fillCaches(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf,
                           InteractionVolume& iv,
                           FluxVarsCacheContainer& fluxVarsCacheContainer)
    {
        // lambda function to get the permeability tensor
        auto getK = [&problem](const Element& element,
                               const VolumeVariables& volVars,
                               const SubControlVolume& scv)
        { return volVars.permeability(); };

        // solve the local system subject to the permeability
        iv.solveLocalSystem(getK);

        // update the transmissibilities etc for each phase
        const auto scvfIdx = scvf.index();
        const auto scvfLocalFaceData = iv.getLocalFaceData(scvf);
        auto& scvfCache = fluxVarsCacheContainer[scvfIdx];
        scvfCache.updateAdvection(iv, scvf, scvfLocalFaceData);

        //! Update the transmissibilities in the other scvfs of the interaction volume
        for (auto scvfIdxJ : iv.globalScvfs())
        {
            if (scvfIdxJ != scvfIdx)
            {
                // update cache of scvfJ
                const auto& scvfJ = fvGeometry.scvf(scvfIdxJ);
                const auto scvfJLocalFaceData = iv.getLocalFaceData(scvfJ);
                fluxVarsCacheContainer[scvfIdxJ].updateAdvection(iv, scvfJ, scvfJLocalFaceData);
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
    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
    static constexpr bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);

public:
    //! function to fill the flux var caches
    template<class FluxVarsCacheContainer, class InteractionVolume>
    static void fillCaches(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf,
                           InteractionVolume& iv,
                           FluxVarsCacheContainer& fluxVarsCacheContainer)
    {
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
                const auto scvfLocalFaceData = iv.getLocalFaceData(scvf);
                auto& scvfCache = fluxVarsCacheContainer[scvfIdx];
                scvfCache.updateDiffusion(iv, scvf, scvfLocalFaceData, phaseIdx, compIdx);

                //! Update the caches in the other scvfs of the interaction volume
                for (auto scvfIdxJ : iv.globalScvfs())
                {
                    if (scvfIdxJ != scvfIdx)
                    {
                        // store scvf pointer and local face data
                        const auto& scvfJ = fvGeometry.scvf(scvfIdxJ);
                        const auto scvfJLocalFaceData = iv.getLocalFaceData(scvfJ);

                        // update cache
                        fluxVarsCacheContainer[scvfIdxJ].updateDiffusion(iv, scvfJ, scvfJLocalFaceData, phaseIdx, compIdx);
                    }
                }
            }
        }
    }
};

//! Fills the quantities for heat conduction in the caches
template<class TypeTag>
class CCMpfaHeatConductionCacheFiller
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);

    using Element = typename GridView::template Codim<0>::Entity;

    static const int energyEqIdx = Indices::energyEqIdx;
    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static const bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);

public:
    //! function to fill the flux var caches
    template<class FluxVarsCacheContainer, class InteractionVolume>
    static void fillCaches(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf,
                           InteractionVolume& iv,
                           FluxVarsCacheContainer& fluxVarsCacheContainer)
    {
        // lambda function to get the thermal conductivity
        auto getLambda = [&problem, &fvGeometry](const Element& element,
                                                 const VolumeVariables& volVars,
                                                 const SubControlVolume& scv)
        { return ThermalConductivityModel::effectiveThermalConductivity(volVars,
                                                                        problem.spatialParams(),
                                                                        element,
                                                                        fvGeometry,
                                                                        scv); };

        // solve for the transmissibilities
        iv.solveLocalSystem(getLambda);

        // update the caches
        const auto scvfIdx = scvf.index();
        const auto scvfLocalFaceData = iv.getLocalFaceData(scvf);
        auto& scvfCache = fluxVarsCacheContainer[scvfIdx];
        scvfCache.updateHeatConduction(iv, scvf, scvfLocalFaceData);

        //! Update the caches in the other scvfs of the interaction volume
        for (auto scvfIdxJ : iv.globalScvfs())
        {
            if (scvfIdxJ != scvfIdx)
            {
                // store scvf pointer and local face data
                const auto& scvfJ = fvGeometry.scvf(scvfIdxJ);
                const auto scvfJLocalFaceData = iv.getLocalFaceData(scvfJ);

                // update cache
                fluxVarsCacheContainer[scvfIdxJ].updateHeatConduction(iv, scvfJ, scvfJLocalFaceData);
            }
        }
    }
};



} // end namespace

#endif
