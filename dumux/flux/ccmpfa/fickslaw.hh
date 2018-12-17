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
 * \ingroup CCMpfaDiscretization
 * \brief Fick's law for cell-centered finite volume schemes with multi-point flux approximation
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FICKS_LAW_HH

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>

namespace Dumux
{
//! forward declaration of the method-specific implemetation
template<class TypeTag, DiscretizationMethod discMethod>
class FicksLawImplementation;

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Fick's law for cell-centered finite volume schemes with multi-point flux approximation
 */
template <class TypeTag>
class FicksLawImplementation<TypeTag, DiscretizationMethod::ccmpfa>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
    using BalanceEqOpts = GetPropType<TypeTag, Properties::BalanceEqOpts>;

    static constexpr int numComponents = GetPropType<TypeTag, Properties::ModelTraits>::numFluidComponents();
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
            if (fvGeometry.fvGridGeometry().vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
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
        static constexpr bool considerSecondaryIVs = FVGridGeometry::MpfaHelper::considerSecondaryIVs();
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

    // state the type for the corresponding cache and its filler
    using Cache = MpfaFicksLawCache;

    //! Compute the diffusive flux across an scvf
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
            if(compIdx == FluidSystem::getMainComponent(phaseIdx))
                continue;

            // get the scaling factor for the effective diffusive fluxes
            const auto effFactor = computeEffectivityFactor(elemVolVars, scvf, phaseIdx);

            // if factor is zero, the flux will end up zero anyway
            if (effFactor == 0.0)
                continue;

            // calculate the density at the interface
            const auto rho = interpolateDensity(elemVolVars, scvf, phaseIdx);

            // compute the flux
            if (fluxVarsCache.usesSecondaryIv())
                componentFlux[compIdx] = rho*effFactor*computeVolumeFlux(problem,
                                                                         fluxVarsCache,
                                                                         fluxVarsCache.diffusionSecondaryDataHandle(),
                                                                         phaseIdx, compIdx);
            else
                componentFlux[compIdx] = rho*effFactor*computeVolumeFlux(problem,
                                                                         fluxVarsCache,
                                                                         fluxVarsCache.diffusionPrimaryDataHandle(),
                                                                         phaseIdx, compIdx);
        }

        // accumulate the phase component flux
        for(int compIdx = 0; compIdx < numComponents; compIdx++)
            if(compIdx != FluidSystem::getMainComponent(phaseIdx) && BalanceEqOpts::mainComponentIsBalanced(phaseIdx) && !FluidSystem::isTracerFluidSystem())
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
            Scalar rho = elemVolVars[scvf.insideScvIdx()].molarDensity(phaseIdx);
            for (const auto outsideIdx : scvf.outsideScvIndices())
                rho += elemVolVars[outsideIdx].molarDensity(phaseIdx);
            return rho/(scvf.outsideScvIndices().size()+1);
        }
        else
            return elemVolVars[scvf.outsideScvIdx()].molarDensity(phaseIdx);
    }

    //! Here we want to calculate the factors with which the diffusion coefficient has to be
    //! scaled to get the effective diffusivity. For this we use the effective diffusivity with
    //! a diffusion coefficient of 1.0 as input. Then we scale the transmissibilites during flux
    //! calculation (above) with the harmonic average of the two factors
    static Scalar computeEffectivityFactor(const ElementVolumeVariables& elemVolVars,
                                           const SubControlVolumeFace& scvf,
                                           const unsigned int phaseIdx)
    {
        using EffDiffModel = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;

        // use the harmonic mean between inside and outside
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto factor = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(),
                                                               insideVolVars.saturation(phaseIdx),
                                                               /*Diffusion coefficient*/ 1.0);

        if (!scvf.boundary())
        {
            // interpret outside factor as arithmetic mean
            Scalar outsideFactor = 0.0;
            for (const auto outsideIdx : scvf.outsideScvIndices())
            {
                const auto& outsideVolVars = elemVolVars[outsideIdx];
                outsideFactor += EffDiffModel::effectiveDiffusivity(outsideVolVars.porosity(),
                                                                    outsideVolVars.saturation(phaseIdx),
                                                                    /*Diffusion coefficient*/ 1.0);
            }
            outsideFactor /= scvf.outsideScvIndices().size();

            // use the harmonic mean of the two
            return harmonicMean(factor, outsideFactor);
        }

        return factor;
    }
};

} // end namespace

#endif
