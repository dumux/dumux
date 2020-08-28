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
 * \ingroup CCMpfaFlux
 * \brief Fick's law for cell-centered finite volume schemes with multi-point flux approximation
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FICKS_LAW_HH

#include <dune/common/fvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>

#include <dumux/flux/fickiandiffusioncoefficients.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

//! forward declaration of the method-specific implemetation
template<class TypeTag, DiscretizationMethod discMethod, ReferenceSystemFormulation referenceSystem>
class FicksLawImplementation;

/*!
 * \ingroup CCMpfaFlux
 * \brief Fick's law for cell-centered finite volume schemes with multi-point flux approximation
 */
template <class TypeTag, ReferenceSystemFormulation referenceSystem>
class FicksLawImplementation<TypeTag, DiscretizationMethod::ccmpfa, referenceSystem>
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
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
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
            // get interaction volume related data from the filler class & upate the cache
            if (fvGeometry.gridGeometry().vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
                scvfFluxVarsCache.updateDiffusion(fluxVarsCacheFiller.secondaryInteractionVolume(),
                                                  fluxVarsCacheFiller.secondaryIvLocalFaceData(),
                                                  fluxVarsCacheFiller.secondaryIvDataHandle(),
                                                  phaseIdx, compIdx);
            else
                scvfFluxVarsCache.updateDiffusion(fluxVarsCacheFiller.primaryInteractionVolume(),
                                                  fluxVarsCacheFiller.primaryIvLocalFaceData(),
                                                  fluxVarsCacheFiller.primaryIvDataHandle(),
                                                  phaseIdx, compIdx);
        }
    };

    //! The cache used in conjunction with the mpfa Fick's Law
    class MpfaFicksLawCache
    {
        using DualGridNodalIndexSet = GetPropType<TypeTag, Properties::DualGridNodalIndexSet>;
        using Stencil = typename DualGridNodalIndexSet::NodalGridStencilType;

        static constexpr int numPhases = GetPropType<TypeTag, Properties::ModelTraits>::numFluidPhases();
        static constexpr bool considerSecondaryIVs = GridGeometry::MpfaHelper::considerSecondaryIVs();
        using PrimaryDataHandle = typename ElementFluxVariablesCache::PrimaryIvDataHandle::DiffusionHandle;
        using SecondaryDataHandle = typename ElementFluxVariablesCache::SecondaryIvDataHandle::DiffusionHandle;

        //! sets the pointer to the data handle (overload for secondary data handles)
        template< bool doSecondary = considerSecondaryIVs, std::enable_if_t<doSecondary, int> = 0 >
        void setHandlePointer_(const SecondaryDataHandle& dataHandle)
        { secondaryHandlePtr_ = &dataHandle; }

        //! sets the pointer to the data handle (overload for primary data handles)
        void setHandlePointer_(const PrimaryDataHandle& dataHandle)
        { primaryHandlePtr_ = &dataHandle; }

    public:
        // export filler type
        using Filler = MpfaFicksLawCacheFiller;

        /*!
         * \brief Update cached objects (transmissibilities).
         *        This is used for updates with primary interaction volumes.
         *
         * \param iv The interaction volume this scvf is embedded in
         * \param localFaceData iv-local info on this scvf
         * \param dataHandle Transmissibility matrix & gravity data of this iv
         */
        template<class IV, class LocalFaceData, class DataHandle>
        void updateDiffusion(const IV& iv,
                             const LocalFaceData& localFaceData,
                             const DataHandle& dataHandle,
                             unsigned int phaseIdx, unsigned int compIdx)
        {
            stencil_[phaseIdx][compIdx] = &iv.stencil();
            switchFluxSign_[phaseIdx][compIdx] = localFaceData.isOutsideFace();
            setHandlePointer_(dataHandle.diffusionHandle());
        }

        //! The stencils corresponding to the transmissibilities
        const Stencil& diffusionStencil(unsigned int phaseIdx, unsigned int compIdx) const
        { return *stencil_[phaseIdx][compIdx]; }

        //! The corresponding data handles
        const PrimaryDataHandle& diffusionPrimaryDataHandle() const { return *primaryHandlePtr_; }
        const SecondaryDataHandle& diffusionSecondaryDataHandle() const { return *secondaryHandlePtr_; }

        //! Returns whether or not this scvf is an "outside" face in the scope of the iv.
        bool diffusionSwitchFluxSign(unsigned int phaseIdx, unsigned int compIdx) const
        { return switchFluxSign_[phaseIdx][compIdx]; }


    private:
        //! phase-/component- specific data
        std::array< std::array<bool, numComponents>, numPhases > switchFluxSign_;
        std::array< std::array<const Stencil*, numComponents>, numPhases > stencil_;

        //! pointers to the corresponding iv-data handles
        const PrimaryDataHandle* primaryHandlePtr_;
        const SecondaryDataHandle* secondaryHandlePtr_;
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::ccmpfa;

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

            // calculate the density at the interface
            const auto rho = interpolateDensity(elemVolVars, scvf, phaseIdx);

            // compute the flux
            if (fluxVarsCache.usesSecondaryIv())
                componentFlux[compIdx] = rho*computeVolumeFlux(problem,
                                                               fluxVarsCache,
                                                               fluxVarsCache.diffusionSecondaryDataHandle(),
                                                               phaseIdx, compIdx);
            else
                componentFlux[compIdx] = rho*computeVolumeFlux(problem,
                                                               fluxVarsCache,
                                                               fluxVarsCache.diffusionPrimaryDataHandle(),
                                                               phaseIdx, compIdx);
        }

        // accumulate the phase component flux
        for (int compIdx = 0; compIdx < numComponents; compIdx++)
            if constexpr (!FluidSystem::isTracerFluidSystem())
                if (compIdx != FluidSystem::getMainComponent(phaseIdx) && BalanceEqOpts::mainComponentIsBalanced(phaseIdx))
                    componentFlux[FluidSystem::getMainComponent(phaseIdx)] -= componentFlux[compIdx];

        return componentFlux;
    }

private:
    template< class Problem, class FluxVarsCache, class DataHandle >
    static Scalar computeVolumeFlux(const Problem& problem,
                                    const FluxVarsCache& cache,
                                    const DataHandle& dataHandle,
                                    int phaseIdx, int compIdx)
    {
        dataHandle.setPhaseIndex(phaseIdx);
        dataHandle.setComponentIndex(compIdx);

        const bool switchSign = cache.diffusionSwitchFluxSign(phaseIdx, compIdx);

        const auto localFaceIdx = cache.ivLocalFaceIndex();
        const auto idxInOutside = cache.indexInOutsideFaces();
        const auto& xj = dataHandle.uj();
        const auto& tij = dim == dimWorld ? dataHandle.T()[localFaceIdx]
                                          : (!switchSign ? dataHandle.T()[localFaceIdx]
                                                         : dataHandle.tijOutside()[localFaceIdx][idxInOutside]);
        Scalar scvfFlux = tij*xj;

        // switch the sign if necessary
        if (switchSign)
            scvfFlux *= -1.0;

        return scvfFlux;
    }

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
