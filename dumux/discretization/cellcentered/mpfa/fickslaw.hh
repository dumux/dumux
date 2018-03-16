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
#include <dumux/discretization/methods.hh>

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
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using BalanceEqOpts = typename GET_PROP_TYPE(TypeTag, BalanceEqOpts);

    static constexpr int numComponents = GET_PROP_TYPE(TypeTag, ModelTraits)::numComponents();
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
                                                  scvf, phaseIdx, compIdx);
            else
                scvfFluxVarsCache.updateDiffusion(fluxVarsCacheFiller.primaryInteractionVolume(),
                                                  fluxVarsCacheFiller.primaryIvLocalFaceData(),
                                                  fluxVarsCacheFiller.primaryIvDataHandle(),
                                                  scvf, phaseIdx, compIdx);
        }
    };

    //! The cache used in conjunction with the mpfa Fick's Law
    class MpfaFicksLawCache
    {
        using DualGridNodalIndexSet = typename GET_PROP_TYPE(TypeTag, DualGridNodalIndexSet);
        using Stencil = typename DualGridNodalIndexSet::NodalGridStencilType;

        using MpfaHelper = typename FVGridGeometry::MpfaHelper;
        static constexpr bool considerSecondaryIVs = MpfaHelper::considerSecondaryIVs();

        using PrimaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
        using PrimaryIvLocalFaceData = typename PrimaryInteractionVolume::Traits::LocalFaceData;
        using PrimaryIvDataHandle = typename ElementFluxVariablesCache::PrimaryIvDataHandle;
        using PrimaryIvCellVector = typename PrimaryInteractionVolume::Traits::MatVecTraits::CellVector;
        using PrimaryIvTij = typename PrimaryInteractionVolume::Traits::MatVecTraits::TMatrix::row_type;

        using SecondaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, SecondaryInteractionVolume);
        using SecondaryIvLocalFaceData = typename SecondaryInteractionVolume::Traits::LocalFaceData;
        using SecondaryIvDataHandle = typename ElementFluxVariablesCache::SecondaryIvDataHandle;
        using SecondaryIvCellVector = typename SecondaryInteractionVolume::Traits::MatVecTraits::CellVector;
        using SecondaryIvTij = typename SecondaryInteractionVolume::Traits::MatVecTraits::TMatrix::row_type;

        static constexpr int dim = GridView::dimension;
        static constexpr int dimWorld = GridView::dimensionworld;
        static constexpr int numPhases = GET_PROP_TYPE(TypeTag, ModelTraits)::numPhases();

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
         * \param scvf The sub-control volume face
         */
        void updateDiffusion(const PrimaryInteractionVolume& iv,
                             const PrimaryIvLocalFaceData& localFaceData,
                             const PrimaryIvDataHandle& dataHandle,
                             const SubControlVolumeFace &scvf,
                             unsigned int phaseIdx, unsigned int compIdx)
        {
            stencil_[phaseIdx][compIdx] = &iv.stencil();
            switchFluxSign_[phaseIdx][compIdx] = localFaceData.isOutside();

            // store pointer to the mole fraction vector of this iv
            primaryXj_[phaseIdx][compIdx] = &dataHandle.moleFractions(phaseIdx, compIdx);

            const auto ivLocalIdx = localFaceData.ivLocalScvfIndex();
            if (dim == dimWorld)
                primaryTij_[phaseIdx][compIdx] = &dataHandle.diffusionT()[ivLocalIdx];
            else
                primaryTij_[phaseIdx][compIdx] = localFaceData.isOutside() ? &dataHandle.diffusionTout()[ivLocalIdx][localFaceData.scvfLocalOutsideScvfIndex()]
                                                                           : &dataHandle.diffusionT()[ivLocalIdx];
        }

        /*!
         * \brief Update cached objects (transmissibilities).
         *        This is used for updates with secondary interaction volumes.
         *
         * \param iv The interaction volume this scvf is embedded in
         * \param localFaceData iv-local info on this scvf
         * \param dataHandle Transmissibility matrix & gravity data of this iv
         * \param scvf The sub-control volume face
         */
        template< bool doSecondary = considerSecondaryIVs, std::enable_if_t<doSecondary, int > = 0 >
        void updateDiffusion(const SecondaryInteractionVolume& iv,
                             const SecondaryIvLocalFaceData& localFaceData,
                             const SecondaryIvDataHandle& dataHandle,
                             const SubControlVolumeFace &scvf,
                             unsigned int phaseIdx, unsigned int compIdx)
        {
            stencil_[phaseIdx][compIdx] = &iv.stencil();
            switchFluxSign_[phaseIdx][compIdx] = localFaceData.isOutside();

            // store pointer to the mole fraction vector of this iv
            secondaryXj_[phaseIdx][compIdx] = &dataHandle.moleFractions(phaseIdx, compIdx);

            const auto ivLocalIdx = localFaceData.ivLocalScvfIndex();
            if (dim == dimWorld)
                secondaryTij_[phaseIdx][compIdx] = &dataHandle.diffusionT()[ivLocalIdx];
            else
                secondaryTij_[phaseIdx][compIdx] = localFaceData.isOutside() ? &dataHandle.diffusionTout()[ivLocalIdx][localFaceData.scvfLocalOutsideScvfIndex()]
                                                                             : &dataHandle.diffusionT()[ivLocalIdx];
        }

        //! In the interaction volume-local system of eq we have one unknown per face.
        //! On scvfs on this face, but in "outside" (neighbor) elements of it, we have
        //! to take the negative value of the fluxes due to the flipped normal vector.
        //! This function returns whether or not this scvf is an "outside" face in the iv.
        bool diffusionSwitchFluxSign(unsigned int phaseIdx, unsigned int compIdx) const
        { return switchFluxSign_[phaseIdx][compIdx]; }

        //! Coefficients for the cell (& Dirichlet) unknowns in flux expressions (primary type)
        const PrimaryIvTij& diffusionTijPrimaryIv(unsigned int phaseIdx, unsigned int compIdx) const
        { return *primaryTij_[phaseIdx][compIdx]; }

        //! Coefficients for the cell (& Dirichlet) unknowns in flux expressions (secondary type)
        const SecondaryIvTij& diffusionTijSecondaryIv(unsigned int phaseIdx, unsigned int compIdx) const
        { return *secondaryTij_[phaseIdx][compIdx]; }

        //! The cell (& maybe Dirichlet) mole fractions within this interaction volume (primary type)
        const PrimaryIvCellVector& moleFractionsPrimaryIv(unsigned int phaseIdx, unsigned int compIdx) const
        { return *primaryXj_[phaseIdx][compIdx]; }

        //! The cell (& maybe Dirichlet) mole fractions within this interaction volume (secondary type)
        const SecondaryIvCellVector& moleFractionsSecondaryIv(unsigned int phaseIdx, unsigned int compIdx) const
        { return *secondaryXj_[phaseIdx][compIdx]; }

        //! The stencils corresponding to the transmissibilities
        const Stencil& diffusionStencil(unsigned int phaseIdx, unsigned int compIdx) const
        { return *stencil_[phaseIdx][compIdx]; }

    private:
        std::array< std::array<bool, numComponents>, numPhases > switchFluxSign_;

        //! The stencils, i.e. the grid indices j
        std::array< std::array<const Stencil*, numComponents>, numPhases > stencil_;

        //! The transmissibilities such that f = Tij*xj
        std::array< std::array<const PrimaryIvCellVector*, numComponents>, numPhases > primaryTij_;
        std::array< std::array<const SecondaryIvCellVector*, numComponents>, numPhases > secondaryTij_;

        //! The interaction-volume wide mole fractions xj
        std::array< std::array<const PrimaryIvCellVector*, numComponents>, numPhases > primaryXj_;
        std::array< std::array<const SecondaryIvCellVector*, numComponents>, numPhases > secondaryXj_;
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

            // calculate Tij*xj
            Scalar flux;
            if (fluxVarsCache.usesSecondaryIv())
            {
                const auto& tij = fluxVarsCache.diffusionTijSecondaryIv(phaseIdx, compIdx);
                const auto& xj = fluxVarsCache.moleFractionsSecondaryIv(phaseIdx, compIdx);
                flux = tij*xj;
            }
            else
            {
                const auto& tij = fluxVarsCache.diffusionTijPrimaryIv(phaseIdx, compIdx);
                const auto& xj = fluxVarsCache.moleFractionsPrimaryIv(phaseIdx, compIdx);
                flux = tij*xj;
            }

            if (fluxVarsCache.diffusionSwitchFluxSign(phaseIdx, compIdx))
                flux *= -1.0;

            componentFlux[compIdx] = flux*rho*effFactor;
        }

        // accumulate the phase component flux
        for(int compIdx = 0; compIdx < numComponents; compIdx++)
            if(compIdx != FluidSystem::getMainComponent(phaseIdx) && BalanceEqOpts::mainComponentIsBalanced(phaseIdx) && !FluidSystem::isTracerFluidSystem())
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
