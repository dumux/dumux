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
 * \brief This file contains the data which is required to calculate
 *        molar and mass fluxes of a component in a fluid phase over a face of a finite volume by means
 *        of Fick's Law for cell-centered MPFA models.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FICKS_LAW_HH

#include <memory>

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

#include <dumux/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup CCMpfaFicksLaw
 * \brief Specialization of Fick's Law for the CCMpfa method.
 */
template <class TypeTag>
class FicksLawImplementation<TypeTag, DiscretizationMethods::CCMpfa>
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

    // Always use the dynamic type for vectors (compatibility with the boundary)
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using CoefficientVector = typename BoundaryInteractionVolume::Vector;

    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static constexpr int numComponents = GET_PROP_VALUE(TypeTag,NumComponents);
    static constexpr bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static constexpr bool enableInteriorBoundaries = GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries);
    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;

    //! The cache used in conjunction with the mpfa Fick's Law
    class MpfaFicksLawCache
    {
         // We always use the dynamic types here to be compatible on the boundary
        using Stencil = typename BoundaryInteractionVolume::GlobalIndexSet;
        using PositionVector = typename BoundaryInteractionVolume::PositionVector;

    public:
        //! The constructor. Initializes the Neumann flux to zero
        MpfaFicksLawCache() { componentNeumannFluxes_.fill(0.0); }

        // update cached objects for the diffusive fluxes
        template<typename InteractionVolume>
        void updateDiffusion(const InteractionVolume& iv, const SubControlVolumeFace &scvf,
                             unsigned int phaseIdx, unsigned int compIdx)
        {
            const auto& localFaceData = iv.getLocalFaceData(scvf);
            diffusionTij_[phaseIdx][compIdx] = iv.getTransmissibilities(localFaceData);
            // copy the stencil only for the first call
            if (phaseIdx == 0 && compIdx == 1)
                diffusionVolVarsStencil_ = iv.volVarsStencil();

            //! For compositional models, we associate neumann fluxes with the phases (main components)
            //! This is done in the AdvectionCache. However, in single-phasic models we solve the phase AND
            //! the component mass balance equations. Thus, in this case we have diffusive neumann contributions.
            //! we assume compIdx = eqIdx
            if (numPhases == 1 && phaseIdx != compIdx)
                componentNeumannFluxes_[compIdx] = iv.getNeumannFlux(localFaceData, compIdx);

        }

        //! Returns the volume variables indices necessary for diffusive flux
        //! computation. This includes all participating boundary volume variables
        //! and it can be different for the phases & components.
        const Stencil& diffusionVolVarsStencil(unsigned int phaseIdx, unsigned int compIdx) const
        { return diffusionVolVarsStencil_; }

        //! Returns the transmissibilities associated with the volume variables
        //! This can be different for the phases & components.
        const CoefficientVector& diffusionTij(unsigned int phaseIdx, unsigned int compIdx) const
        { return diffusionTij_[phaseIdx][compIdx]; }

        //! If the useTpfaBoundary property is set to false, the boundary conditions
        //! are put into the local systems leading to possible contributions on all faces
        Scalar componentNeumannFlux(unsigned int compIdx) const
        {
            assert(numPhases == 1);
            return componentNeumannFluxes_[compIdx];
        }

    private:
        // Quantities associated with molecular diffusion
        Stencil diffusionVolVarsStencil_;
        std::array< std::array<CoefficientVector, numComponents>, numPhases> diffusionTij_;

        // diffusive neumann flux for single-phasic models
        std::array<Scalar, numComponents> componentNeumannFluxes_;
    };

    //! Class that fills the cache corresponding to mpfa Darcy's Law
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
            if (problem.model().globalFvGeometry().isInBoundaryInteractionVolume(scvf))
                scvfFluxVarsCache.updateDiffusion(fluxVarsCacheFiller.boundaryInteractionVolume(), scvf, phaseIdx, compIdx);
            else
                scvfFluxVarsCache.updateDiffusion(fluxVarsCacheFiller.interactionVolume(), scvf, phaseIdx, compIdx);
        }
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCMpfa;

    // state the type for the corresponding cache and its filler
    using Cache = MpfaFicksLawCache;
    using CacheFiller = MpfaFicksLawCacheFiller;

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
            if(compIdx == phaseIdx)
              continue;

            const auto& fluxVarsCache = elemFluxVarsCache[scvf];
            const auto& volVarsStencil = fluxVarsCache.diffusionVolVarsStencil(phaseIdx, compIdx);
            const auto& tij = fluxVarsCache.diffusionTij(phaseIdx, compIdx);

            const bool isInteriorBoundary = enableInteriorBoundaries && fluxVarsCache.isInteriorBoundary();
            // For interior Neumann boundaries when using Tpfa for Neumann boundary conditions, we simply
            // return the user-specified flux. Note that for compositional models we attribute the influxes
            // to the major components, thus we do it per phase in Darcy's law. However, for single-phasic models
            // wesolve the phase mass balance equation AND the transport equation, thus, in that case we incorporate
            // the Neumann BCs here. We assume compIdx = eqIdx.
            // Note that this way of including interior Neumann fluxes fails for mpnc models where n != m.
            if (numPhases == 1
                && isInteriorBoundary
                && useTpfaBoundary
                && fluxVarsCache.interiorBoundaryDataSelf().faceType() == MpfaFaceTypes::interiorNeumann)
                componentFlux[compIdx] = scvf.area()*
                                         elemVolVars[scvf.insideScvIdx()].extrusionFactor()*
                                         problem.neumann(element, fvGeometry, elemVolVars, scvf)[compIdx];


            // get the scaling factor for the effective diffusive fluxes
            const auto effFactor = Implementation::computeEffectivityFactor(fvGeometry, elemVolVars, scvf, fluxVarsCache, phaseIdx, isInteriorBoundary);

            // if factor is zero, the flux will end up zero anyway
            if (effFactor == 0.0)
            {
                componentFlux[compIdx] = 0.0;
                continue;
            }

            // lambda functions depending on if we use mole or mass fractions
            auto getX = [phaseIdx, compIdx] (const auto& volVars)
            { return volVars.moleFraction(phaseIdx, compIdx); };

            auto getRho = [phaseIdx] (const auto& volVars)
            { return volVars.molarDensity(phaseIdx); };

            // calculate the density at the interface
            const auto rho = Implementation::interpolateDensity(fvGeometry, elemVolVars, scvf, fluxVarsCache, getRho, isInteriorBoundary);

            // calculate Tij*xj
            Scalar flux(0.0);
            unsigned int localIdx = 0;
            for (const auto volVarIdx : volVarsStencil)
                flux += tij[localIdx++]*getX(elemVolVars[volVarIdx]);

            // if no interior boundaries are present, return effective mass flux
            if (!enableInteriorBoundaries)
                componentFlux[compIdx] = useTpfaBoundary ? flux*rho*effFactor : flux*rho*effFactor + fluxVarsCache.componentNeumannFlux(compIdx);

            else
            {
              // Handle interior boundaries
              flux += Implementation::computeInteriorBoundaryContribution(fvGeometry, elemVolVars, fluxVarsCache, getX, phaseIdx, compIdx);

              // return overall resulting flux
              componentFlux[compIdx] = useTpfaBoundary ? flux*rho*effFactor : flux*rho*effFactor + fluxVarsCache.componentNeumannFlux(compIdx);
            }
        }

        // accumulate the phase component flux
        for(int compIdx = 0; compIdx < numComponents; compIdx++)
          if(compIdx != phaseIdx)
            componentFlux[phaseIdx] -= componentFlux[compIdx];

        return componentFlux;
    }

    template<typename GetRhoFunction>
    static Scalar interpolateDensity(const FVElementGeometry& fvGeometry,
                                     const ElementVolumeVariables& elemVolVars,
                                     const SubControlVolumeFace& scvf,
                                     const FluxVariablesCache& fluxVarsCache,
                                     const GetRhoFunction& getRho,
                                     const bool isInteriorBoundary)
    {

        // maybe use the density of the interior BC on the facet
        if (isInteriorBoundary)
        {
            const auto& data = fluxVarsCache.interiorBoundaryDataSelf();
            if (data.faceType() == MpfaFaceTypes::interiorDirichlet)
                return getRho(data.facetVolVars(fvGeometry));
        }

        // use arithmetic mean of the densities around the scvf
        if (!scvf.boundary())
        {
            Scalar rho = getRho(elemVolVars[scvf.insideScvIdx()]);
            for (auto outsideIdx : scvf.outsideScvIndices())
                rho += getRho(elemVolVars[outsideIdx]);
            return rho/(scvf.outsideScvIndices().size()+1);
        }
        else
            return getRho(elemVolVars[scvf.outsideScvIdx()]);
    }

    //! Here we want to calculate the factors with which the diffusion coefficient has to be
    //! scaled to get the effective diffusivity. For this we use the effective diffusivity with
    //! a diffusion coefficient of 1.0 as input. Then we scale the transmissibilites during flux
    //! calculation (above) with the harmonic average of the two factors
    static Scalar computeEffectivityFactor(const FVElementGeometry& fvGeometry,
                                           const ElementVolumeVariables& elemVolVars,
                                           const SubControlVolumeFace& scvf,
                                           const FluxVariablesCache& fluxVarsCache,
                                           const unsigned int phaseIdx,
                                           const bool isInteriorBoundary)
    {
        // Treat interior boundaries differently
        if (isInteriorBoundary)
        {
            const auto& data = fluxVarsCache.interiorBoundaryDataSelf();
            // use harmonic mean between the interior and the facet volvars
            if (data.faceType() == MpfaFaceTypes::interiorDirichlet)
            {
                const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
                const auto factor = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(),
                                                                       insideVolVars.saturation(phaseIdx),
                                                                       /*Diffusion coefficient*/ 1.0);

                const auto facetVolVars = data.facetVolVars(fvGeometry);
                const auto outsideFactor = EffDiffModel::effectiveDiffusivity(facetVolVars.porosity(),
                                                                              facetVolVars.saturation(phaseIdx),
                                                                              /*Diffusion coefficient*/ 1.0);

                // check if we divide by zero
                if (factor*outsideFactor <= 0.0)
                    return 0.0;
                return harmonicMean(factor, outsideFactor);
            }
        }

        // use the harmonic mean between inside and outside
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto factor = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(),
                                                               insideVolVars.saturation(phaseIdx),
                                                               /*Diffusion coefficient*/ 1.0);

        if (!scvf.boundary())
        {
            // interpret outside factor as arithmetic mean
            Scalar outsideFactor = 0.0;
            for (auto outsideIdx : scvf.outsideScvIndices())
            {
                const auto& outsideVolVars = elemVolVars[outsideIdx];
                outsideFactor += EffDiffModel::effectiveDiffusivity(outsideVolVars.porosity(),
                                                                    outsideVolVars.saturation(phaseIdx),
                                                                    /*Diffusion coefficient*/ 1.0);
            }
            outsideFactor /= scvf.outsideScvIndices().size();

            // check if we divide by zero
            if (factor*outsideFactor <= 0.0)
                return 0.0;
            return harmonicMean(factor, outsideFactor);
        }

        return factor;
    }

    template<typename GetXFunction>
    static Scalar computeInteriorBoundaryContribution(const FVElementGeometry& fvGeometry,
                                                      const ElementVolumeVariables& elemVolVars,
                                                      const FluxVariablesCache& fluxVarsCache,
                                                      const GetXFunction& getX,
                                                      unsigned int phaseIdx, unsigned int compIdx)
    {
        // obtain the transmissibilites associated with all pressures
        const auto& tij = fluxVarsCache.diffusionTij(phaseIdx, compIdx);

        // the interior dirichlet boundaries local indices start after
        // the cell and the domain Dirichlet boundary pressures
        const auto startIdx = fluxVarsCache.diffusionVolVarsStencil(phaseIdx, compIdx).size();

        // add interior Dirichlet boundary contributions
        Scalar flux = 0.0;
        for (auto&& data : fluxVarsCache.interiorBoundaryData())
            if (data.faceType() == MpfaFaceTypes::interiorDirichlet)
                flux += tij[startIdx + data.localIndexInInteractionVolume()]*getX(data.facetVolVars(fvGeometry));

        return flux;
    }
};

} // end namespace

#endif
