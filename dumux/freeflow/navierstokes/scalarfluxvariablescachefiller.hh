// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \brief A helper class to fill the flux variables cache
 */
#ifndef DUMUX_FREEFLOW_SCALAR_FLUXVARIABLESCACHE_FILLER_HH
#define DUMUX_FREEFLOW_SCALAR_FLUXVARIABLESCACHE_FILLER_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/referencesystemformulation.hh>
#include <dumux/common/typetraits/problem.hh>

namespace Dumux {

// forward declaration
template<class Problem, class ModelTraits, bool diffusionIsSolDependent, bool heatConductionIsSolDependent, class DiscretizationMethod>
class FreeFlowScalarFluxVariablesCacheFillerImplementation;

/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables cache filler class for free flow
 *
 * Helps filling the flux variables cache depending several policies
 */
template<class Problem, class ModelTraits, bool diffusionIsSolDependent, bool heatConductionIsSolDependent>
using FreeFlowScalarFluxVariablesCacheFiller = FreeFlowScalarFluxVariablesCacheFillerImplementation<Problem, ModelTraits, diffusionIsSolDependent, heatConductionIsSolDependent, typename ProblemTraits<Problem>::GridGeometry::DiscretizationMethod>;

//! Specialization of the flux variables cache filler for the cell centered tpfa method
template<class Problem, class ModelTraits, bool diffusionIsSolDependent, bool heatConductionIsSolDependent>
class FreeFlowScalarFluxVariablesCacheFillerImplementation<Problem, ModelTraits, diffusionIsSolDependent, heatConductionIsSolDependent, DiscretizationMethods::CCTpfa>
{
    using GridGeometry = typename ProblemTraits<Problem>::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr bool diffusionEnabled = ModelTraits::enableMolecularDiffusion();
    static constexpr bool heatConductionEnabled = ModelTraits::enableEnergyBalance();

public:
    static constexpr bool isSolDependent = (diffusionEnabled && diffusionIsSolDependent) ||
                                           (heatConductionEnabled && heatConductionIsSolDependent);

    //! The constructor. Sets the problem pointer
    FreeFlowScalarFluxVariablesCacheFillerImplementation(const Problem& problem)
    : problemPtr_(&problem) {}

    /*!
     * \brief function to fill the flux variables caches
     *
     * \param fluxVarsCacheContainer Either the element or global flux variables cache
     * \param scvfFluxVarsCache The flux var cache to be updated corresponding to the given scvf
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param elemVolVars The element volume variables
     * \param scvf The corresponding sub-control volume face
     * \param forceUpdateAll if true, forces all caches to be updated (even the solution-independent ones)
     */
    template<class FluxVariablesCacheContainer, class FluxVariablesCache, class ElementVolumeVariables>
    void fill(FluxVariablesCacheContainer& fluxVarsCacheContainer,
              FluxVariablesCache& scvfFluxVarsCache,
              const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars,
              const SubControlVolumeFace& scvf,
              const bool forceUpdateAll = false)
    {
        // fill the physics-related quantities of the caches
        if (forceUpdateAll)
        {
            if constexpr (diffusionEnabled)
                fillDiffusion_(scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf);
            if constexpr (heatConductionEnabled)
                fillHeatConduction_(scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf);
        }
        else
        {
            if constexpr (diffusionEnabled && diffusionIsSolDependent)
                fillDiffusion_(scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf);
            if constexpr (heatConductionEnabled && heatConductionIsSolDependent)
                fillHeatConduction_(scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf);
        }
    }

private:

    const Problem& problem() const
    { return *problemPtr_; }


    //! method to fill the diffusive quantities
    template<class FluxVariablesCache, class ElementVolumeVariables>
    void fillDiffusion_(FluxVariablesCache& scvfFluxVarsCache,
                        const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf)
    {
        using DiffusionType = typename ElementVolumeVariables::VolumeVariables::MolecularDiffusionType;
        using DiffusionFiller = typename DiffusionType::Cache::Filler;
        using FluidSystem = typename ElementVolumeVariables::VolumeVariables::FluidSystem;

        static constexpr int numPhases = ModelTraits::numFluidPhases();
        static constexpr int numComponents = ModelTraits::numFluidComponents();

        // forward to the filler of the diffusive quantities
        if constexpr (FluidSystem::isTracerFluidSystem())
            for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                for (unsigned int compIdx = 0; compIdx < numComponents; ++compIdx)
                    DiffusionFiller::fill(scvfFluxVarsCache, phaseIdx, compIdx, problem(), element, fvGeometry, elemVolVars, scvf, *this);
        else
            for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                for (unsigned int compIdx = 0; compIdx < numComponents; ++compIdx)
                    if (compIdx != FluidSystem::getMainComponent(phaseIdx))
                        DiffusionFiller::fill(scvfFluxVarsCache, phaseIdx, compIdx, problem(), element, fvGeometry, elemVolVars, scvf, *this);
    }

    //! method to fill the quantities related to heat conduction
    template<class FluxVariablesCache, class ElementVolumeVariables>
    void fillHeatConduction_(FluxVariablesCache& scvfFluxVarsCache,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf)
    {
        using HeatConductionType = typename ElementVolumeVariables::VolumeVariables::HeatConductionType;
        using HeatConductionFiller = typename HeatConductionType::Cache::Filler;

        // forward to the filler of the diffusive quantities
        HeatConductionFiller::fill(scvfFluxVarsCache, problem(), element, fvGeometry, elemVolVars, scvf, *this);
    }

    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
