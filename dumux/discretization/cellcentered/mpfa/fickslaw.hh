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
 * \brief This file contains the class which is required to calculate
 *        molar and mass fluxes of a component in a fluid phase over a face of a finite volume by means
 *        of Fick's Law for cell-centered MPFA models.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FICKS_LAW_HH

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, DiscretizationMethods discMethod>
class FicksLawImplementation;

/*!
 * \ingroup Mpfa
 * \brief Specialization of Fick's Law for the CCMpfa method.
 */
template <class TypeTag>
class FicksLawImplementation<TypeTag, DiscretizationMethods::CCMpfa>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using BalanceEqOpts = typename GET_PROP_TYPE(TypeTag, BalanceEqOpts);

    // Always use the dynamic type for vectors (compatibility with the boundary)
    using PrimaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
    using CoefficientVector = typename PrimaryInteractionVolume::Traits::DynamicVector;
    using DataHandle = typename PrimaryInteractionVolume::Traits::DataHandle;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static constexpr int numComponents = GET_PROP_VALUE(TypeTag,NumComponents);
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
          // get interaction volume from the flux vars cache filler & upate the cache
          if (fvGeometry.fvGridGeometry().vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
              scvfFluxVarsCache.updateDiffusion(fluxVarsCacheFiller.secondaryInteractionVolume(),
                                                fluxVarsCacheFiller.dataHandle(),
                                                scvf, phaseIdx, compIdx);
          else
              scvfFluxVarsCache.updateDiffusion(fluxVarsCacheFiller.primaryInteractionVolume(),
                                                fluxVarsCacheFiller.dataHandle(),
                                                scvf, phaseIdx, compIdx);
        }
    };

    //! The cache used in conjunction with the mpfa Fick's Law
    class MpfaFicksLawCache
    {
        // We always use the dynamic types here to be compatible on the boundary
        using Stencil = typename PrimaryInteractionVolume::Traits::DynamicGlobalIndexContainer;
        using DirichletDataContainer = typename PrimaryInteractionVolume::DirichletDataContainer;

    public:
        // export filler type
        using Filler = MpfaFicksLawCacheFiller;

        // update cached objects for the diffusive fluxes
        template<class InteractionVolume>
        void updateDiffusion(const InteractionVolume& iv,
                             const DataHandle& dataHandle,
                             const SubControlVolumeFace &scvf,
                             unsigned int phaseIdx, unsigned int compIdx)
        {
            const auto& localFaceData = iv.getLocalFaceData(scvf);

            // update the quantities that are equal for all phases
            diffusionSwitchFluxSign_[phaseIdx][compIdx] = localFaceData.isOutside();
            diffusionVolVarsStencil_[phaseIdx][compIdx] = &dataHandle.volVarsStencil();
            diffusionDirichletData_[phaseIdx][compIdx] = &dataHandle.dirichletData();

            // the transmissibilities on surface grids have to be obtained from the outside
            if (dim == dimWorld)
                diffusionTij_[phaseIdx][compIdx] = &dataHandle.T()[localFaceData.ivLocalScvfIndex()];
            else
                diffusionTij_[phaseIdx][compIdx] = localFaceData.isOutside() ?
                                                   &dataHandle.outsideTij()[localFaceData.ivLocalOutsideScvfIndex()] :
                                                   &dataHandle.T()[localFaceData.ivLocalScvfIndex()];
        }

        //! Returns the volume variables indices necessary for diffusive flux
        //! computation. This includes all participating boundary volume variables
        //! and it can be different for the phases & components.
        const Stencil& diffusionVolVarsStencil(unsigned int phaseIdx, unsigned int compIdx) const
        { return *diffusionVolVarsStencil_[phaseIdx][compIdx]; }

        //! On faces that are "outside" w.r.t. a face in the interaction volume,
        //! we have to take the negative value of the fluxes, i.e. multiply by -1.0
        bool diffusionSwitchFluxSign(unsigned int phaseIdx, unsigned int compIdx) const
        { return diffusionSwitchFluxSign_[phaseIdx][compIdx]; }

        //! Returns the transmissibilities associated with the volume variables
        //! This can be different for the phases & components.
        const CoefficientVector& diffusionTij(unsigned int phaseIdx, unsigned int compIdx) const
        { return *diffusionTij_[phaseIdx][compIdx]; }

        //! Returns the data on dirichlet boundary conditions affecting
        //! the flux computation on this face
        const DirichletDataContainer& diffusionDirichletData(unsigned int phaseIdx, unsigned int compIdx) const
        { return *diffusionDirichletData_[phaseIdx][compIdx]; }

    private:
        std::array< std::array<bool, numComponents>, numPhases> diffusionSwitchFluxSign_;
        std::array< std::array<const Stencil*, numComponents>, numPhases> diffusionVolVarsStencil_;
        std::array< std::array<const DirichletDataContainer*, numComponents>, numPhases> diffusionDirichletData_;
        std::array< std::array<const CoefficientVector*, numComponents>, numPhases> diffusionTij_;
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCMpfa;

    // state the type for the corresponding cache and its filler
    using Cache = MpfaFicksLawCache;

    static ComponentFluxVector flux (const Problem& problem,
                                     const Element& element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVolumeVariables&  elemVolVars,
                                     const SubControlVolumeFace& scvf,
                                     const int phaseIdx,
                                     const ElementFluxVariablesCache& elemFluxVarsCache)
    {
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

            // prepare computations
            Scalar flux(0.0);
            unsigned int i = 0;
            const auto& fluxVarsCache = elemFluxVarsCache[scvf];
            const auto& tij = fluxVarsCache.diffusionTij(phaseIdx, compIdx);

            // calculate Tij*xj
            for (const auto volVarIdx : fluxVarsCache.diffusionVolVarsStencil(phaseIdx, compIdx))
                flux += tij[i++]*elemVolVars[volVarIdx].moleFraction(phaseIdx, compIdx);

            // add contributions from dirichlet BCs
            for (const auto& d : fluxVarsCache.diffusionDirichletData(phaseIdx, compIdx))
                flux += tij[i++]*elemVolVars[d.volVarIndex()].moleFraction(phaseIdx, compIdx);

            componentFlux[compIdx] = fluxVarsCache.diffusionSwitchFluxSign(phaseIdx, compIdx) ? -1.0*flux*rho*effFactor : flux*rho*effFactor;
        }

        // accumulate the phase component flux
        for(int compIdx = 0; compIdx < numComponents; compIdx++)
            if(compIdx != FluidSystem::getMainComponent(phaseIdx) && BalanceEqOpts::mainComponentIsBalanced(phaseIdx) && !FluidSystem::isTracerFluidSystem())
                componentFlux[FluidSystem::getMainComponent(phaseIdx)] -= componentFlux[compIdx];

        return componentFlux;
    }

private:
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
        using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);

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
