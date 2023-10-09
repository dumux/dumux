// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCMpfaFlux
 * \brief Fick's law for cell-centered finite volume schemes with multi-point flux approximation
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FICKS_LAW_HH

#include <dune/common/fvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>

#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>
#include <dumux/discretization/cellcentered/mpfa/interactionvolume.hh>

#include <dumux/flux/fickiandiffusioncoefficients.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

//! forward declaration of the method-specific implementation
template<class TypeTag, class DiscretizationMethod, ReferenceSystemFormulation referenceSystem>
class FicksLawImplementation;

/*!
 * \ingroup CCMpfaFlux
 * \brief Fick's law for cell-centered finite volume schemes with multi-point flux approximation
 */
template <class TypeTag, ReferenceSystemFormulation referenceSystem>
class FicksLawImplementation<TypeTag, DiscretizationMethods::CCMpfa, referenceSystem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using GridFluxVariablesCache = GetPropType<TypeTag, Properties::GridFluxVariablesCache>;
    using ElementFluxVariablesCache = typename GridFluxVariablesCache::LocalView;
    using FluxVariablesCache = typename GridFluxVariablesCache::FluxVariablesCache;
    using BalanceEqOpts = GetPropType<TypeTag, Properties::BalanceEqOpts>;

    static constexpr int numComponents = GetPropType<TypeTag, Properties::ModelTraits>::numFluidComponents();
    static constexpr int numPhases = GetPropType<TypeTag, Properties::ModelTraits>::numFluidPhases();
    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;

    //! Class that fills the cache corresponding to mpfa Fick's Law
    class MpfaFicksLawCacheFiller
    {
    public:
        //! Function to fill an MpfaFicksLawCache of a given scvf
        //! This interface has to be met by any diffusion-related cache filler class
        template<class FluxVariablesCacheFiller>
        static void fill(FluxVariablesCache& scvfFluxVarsCache,
                         unsigned int phaseIdx, unsigned int compIdx,
                         const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace& scvf,
                         const FluxVariablesCacheFiller& fluxVarsCacheFiller)
        {
            // get interaction volume related data from the filler class & update the cache
            scvfFluxVarsCache.updateDiffusion(
                fluxVarsCacheFiller.fluxesPtr(),
                phaseIdx, compIdx,
                fluxVarsCacheFiller.diffusionFluxIds()[phaseIdx][compIdx]
            );
        }
    };

    //! The cache used in conjunction with the mpfa Fick's Law
    class MpfaFicksLawCache
    {
        using Fluxes = CCMpfaFluxes<GridGeometry, Scalar>;
        using FluxId = typename Fluxes::FluxId;

        static constexpr int numPhases = GetPropType<TypeTag, Properties::ModelTraits>::numFluidPhases();
        static constexpr int numComponents = GetPropType<TypeTag, Properties::ModelTraits>::numFluidComponents();

    public:
        // export filler type
        using Filler = MpfaFicksLawCacheFiller;

        /*!
         * \brief Update cached objects (i.e. ptrs to flux computation instance).
         * \param fluxes The flux computation instance
         * \param phaseIdx The index of the phase to update
         * \param compIdx The index of the component to update
         * \param ids Flux id for the diffusive flux of the given phase/component
         */
        void updateDiffusion(const Fluxes* fluxesPtr,
                             unsigned int phaseIdx,
                             unsigned int compIdx,
                             FluxId id)
        {
            fluxes_ = fluxesPtr;
            diffusionFluxIds_[phaseIdx][compIdx] = std::move(id);
        }

        const Fluxes& diffusionFluxes() const
        { return *fluxes_; }

        FluxId diffusionId(int phaseIdx, int compIdx) const
        { return diffusionFluxIds_[phaseIdx][compIdx]; }

    private:
        const Fluxes* fluxes_;
        std::array<std::array<FluxId, numComponents>, numPhases> diffusionFluxIds_;
    };

public:
    using DiscretizationMethod = DiscretizationMethods::CCMpfa;
    // state the discretization method this implementation belongs to
    static constexpr DiscretizationMethod discMethod{};

    //return the reference system
    static constexpr ReferenceSystemFormulation referenceSystemFormulation()
    { return referenceSystem; }

    // state the type for the corresponding cache and its filler
    using Cache = MpfaFicksLawCache;

    // export the diffusion container
    using DiffusionCoefficientsContainer = FickianDiffusionCoefficients<Scalar, numPhases, numComponents>;

    /*!
     * \brief Returns the diffusive fluxes of all components within
     *        a fluid phase across the given sub-control volume face.
     *        The computed fluxes are given in mole/s or kg/s, depending
     *        on the template parameter ReferenceSystemFormulation.
     */
    static ComponentFluxVector flux(const Problem& problem,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables&  elemVolVars,
                                    const SubControlVolumeFace& scvf,
                                    const int phaseIdx,
                                    const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        // obtain this scvf's cache
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        ComponentFluxVector componentFlux(0.0);
        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            if constexpr (!FluidSystem::isTracerFluidSystem())
                if (compIdx == FluidSystem::getMainComponent(phaseIdx))
                    continue;

            const auto rho = interpolateDensity(elemVolVars, scvf, phaseIdx);
            componentFlux[compIdx] = rho*fluxVarsCache.diffusionFluxes().computeFluxFor(
                fluxVarsCache.diffusionId(phaseIdx, compIdx),
                scvf
            );
        }

        // accumulate the phase component flux
        for (int compIdx = 0; compIdx < numComponents; compIdx++)
            if constexpr (!FluidSystem::isTracerFluidSystem())
                if (compIdx != FluidSystem::getMainComponent(phaseIdx) && BalanceEqOpts::mainComponentIsBalanced(phaseIdx))
                    componentFlux[FluidSystem::getMainComponent(phaseIdx)] -= componentFlux[compIdx];

        return componentFlux;
    }

private:
    //! compute the density at branching facets for network grids as arithmetic mean
    static Scalar interpolateDensity(const ElementVolumeVariables& elemVolVars,
                                     const SubControlVolumeFace& scvf,
                                     const unsigned int phaseIdx)
    {
        // use arithmetic mean of the densities around the scvf
        if (!scvf.boundary())
        {
            const Scalar rhoInside = massOrMolarDensity(elemVolVars[scvf.insideScvIdx()], referenceSystem, phaseIdx);

            Scalar rho = rhoInside;
            for (const auto outsideIdx : scvf.outsideScvIndices())
            {
                const Scalar rhoOutside = massOrMolarDensity(elemVolVars[outsideIdx], referenceSystem, phaseIdx);
                rho += rhoOutside;
            }
            return rho/(scvf.outsideScvIndices().size()+1);
        }
        else
            return massOrMolarDensity(elemVolVars[scvf.outsideScvIdx()], referenceSystem, phaseIdx);
    }
};

} // end namespace

#endif
