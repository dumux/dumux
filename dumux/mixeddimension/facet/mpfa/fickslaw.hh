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
 *        of Fick's Law for cell-centered MPFA models in the presence of lower dimensional (coupled) elements
 *        living on the this domain's element facets.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FACET_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FACET_FICKS_LAW_HH

#include <dumux/discretization/cellcentered/mpfa/fickslaw.hh>
#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>
#include <dumux/discretization/cellcentered/mpfa/tensorlambdafactory.hh>

namespace Dumux
{
/*!
 * \ingroup CCMpfaFicksLaw
 * \brief Specialization of Fick's Law for the CCMpfa method with lower dimensional
 *        elements living on the bulk elements' facets.
 */
template <class TypeTag>
class CCMpfaFacetCouplingFicksLaw : public FicksLawImplementation<TypeTag, DiscretizationMethods::CCMpfa>
{
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
    static constexpr bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);

    //! The cache used in conjunction with the mpfa Fick's Law
    class MpfaFacetCouplingFicksLawCache
    {
        static const int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);

        // We always use the dynamic types here to be compatible on the boundary
        using Stencil = typename BoundaryInteractionVolume::GlobalIndexSet;
        using PositionVector = typename BoundaryInteractionVolume::PositionVector;

    public:
        //! The constructor. Initializes the Neumann flux to zero
        MpfaFacetCouplingFicksLawCache() { componentNeumannFluxes_.fill(0.0); }

        // update cached objects for the diffusive fluxes
        template<typename InteractionVolume>
        void updateDiffusion(const InteractionVolume& iv, const SubControlVolumeFace &scvf,
                             unsigned int phaseIdx, unsigned int compIdx)
        {
            const auto& localFaceData = iv.getLocalFaceData(scvf);
            diffusionTij_[phaseIdx][compIdx] = iv.getTransmissibilities(localFaceData);
            // get the stencil only for the first call
            if (phaseIdx == 0 && compIdx == 1)
                diffusionVolVarsStencil_ = iv.volVarsStencil();

            // we will need the neumann flux transformation on interior Neumann boundaries
            diffusionCij_[phaseIdx][compIdx] = iv.getNeumannFluxTransformationCoefficients(localFaceData);

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

        //! Returns the vector of coefficients with which the vector of neumann boundary conditions
        //! has to be multiplied in order to transform them on the scvf this cache belongs to
        const CoefficientVector& diffusionCij(unsigned int phaseIdx, unsigned int compIdx) const
        { return diffusionCij_[phaseIdx][compIdx]; }

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
        std::array< std::array<CoefficientVector, numComponents>, numPhases> diffusionCij_;

        // diffusive neumann flux for single-phasic models
        std::array<Scalar, numComponents> componentNeumannFluxes_;
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCMpfa;

    // state the new type for the corresponding cache
    using Cache = MpfaFacetCouplingFicksLawCache;

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       int phaseIdx, int compIdx,
                       const ElementFluxVariablesCache& elemFluxVarsCache,
                       bool useMoles = true)
    {
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto& volVarsStencil = fluxVarsCache.diffusionVolVarsStencil(phaseIdx, compIdx);
        const auto& tij = fluxVarsCache.diffusionTij(phaseIdx, compIdx);

        const bool isInteriorBoundary = fluxVarsCache.isInteriorBoundary();

        // get the scaling factor for the effective diffusive fluxes
        const auto effFactor = numPhases == 1 ? 1.0 :
                               computeEffectivityFactor(fvGeometry, elemVolVars, scvf, fluxVarsCache, phaseIdx, isInteriorBoundary);

        // if factor is zero, the flux will end up zero anyway
        if (effFactor == 0.0)
            return 0.0;

        // lambda functions depending on if we use mole or mass fractions
        auto getX = [useMoles, phaseIdx, compIdx] (const auto& volVars)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        auto getRho = [useMoles, phaseIdx] (const auto& volVars)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        // calculate the density at the interface
        const auto rho = interpolateDensity(problem, element, fvGeometry, elemVolVars, scvf, fluxVarsCache, getRho, isInteriorBoundary);

        // calculate Tij*xj
        Scalar flux(0.0);
        unsigned int localIdx = 0;
        for (const auto volVarIdx : volVarsStencil)
            flux += tij[localIdx++]*getX(elemVolVars[volVarIdx]);

        // Handle interior boundaries
        flux += computeInteriorBoundaryContribution(fvGeometry, elemVolVars, fluxVarsCache, getX, phaseIdx, compIdx);

        // return overall resulting flux
        return useTpfaBoundary ? flux*rho*effFactor : flux*rho*effFactor + fluxVarsCache.componentNeumannFlux(compIdx);
    }

    template<typename GetRhoFunction>
    static Scalar interpolateDensity(const Problem& problem,
                                     const Element& element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVolumeVariables& elemVolVars,
                                     const SubControlVolumeFace& scvf,
                                     const FluxVariablesCache& fluxVarsCache,
                                     const GetRhoFunction& getRho,
                                     const bool isInteriorBoundary)
    {

        // maybe use the density of the interior BC on the facet
        if (isInteriorBoundary)
            return getRho(fluxVarsCache.interiorBoundaryDataSelf().facetVolVars(fvGeometry));

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
    //! scaled to get the effective diffusivity (for numPhases > 1). For this we use the effective
    //! diffusivity with a diffusion coefficient of 1.0 as input. Then we scale the transmissibilites
    //! during flux calculation (above) with the harmonic average of the two factors. We need to do this
    //! as the mpfa cannot handle diffusion coefficients of zero in the local equation systems.
    static Scalar computeEffectivityFactor(const FVElementGeometry& fvGeometry,
                                           const ElementVolumeVariables& elemVolVars,
                                           const SubControlVolumeFace& scvf,
                                           const FluxVariablesCache& fluxVarsCache,
                                           const unsigned int phaseIdx,
                                           const bool isInteriorBoundary)
    {
        if (isInteriorBoundary)
        {
            // use harmonic mean between the interior and the facet volvars
            const auto& data = fluxVarsCache.interiorBoundaryDataSelf();

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

        // use the harmonic mean between inside and outside
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto factor = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(),
                                                               insideVolVars.saturation(phaseIdx),
                                                               /*Diffusion coefficient*/ 1.0);

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

        // The vector of interior neumann fluxes
        const auto& cij = fluxVarsCache.diffusionCij(phaseIdx, compIdx);
        CoefficientVector facetCouplingFluxes(cij.size(), 0.0);

        // add interior Dirichlet boundary contributions
        Scalar flux = 0.0;
        for (auto&& data : fluxVarsCache.interiorBoundaryData())
        {
            // Add additional Dirichlet fluxes for interior Dirichlet faces
            if (data.faceType() == MpfaFaceTypes::interiorDirichlet)
                flux += tij[startIdx + data.localIndexInInteractionVolume()]*getX(data.facetVolVars(fvGeometry));

            // add neumann contributions
            if (data.faceType() == MpfaFaceTypes::interiorNeumann)
            {
                // get the scvf corresponding to actual interior neumann face
                const auto& curScvf = fvGeometry.scvf(data.scvfIndex());

                if (numPhases > 1)
                {
                    // get the volume variables of the actual interior neumann face
                    const auto facetVolVars = data.facetVolVars(fvGeometry);

                    // calculate "leakage factor"
                    const auto n = curScvf.unitOuterNormal();
                    const auto v = [&] ()
                                    {
                                        auto res = n;
                                        res *= -0.5*facetVolVars.extrusionFactor();
                                        res /= res.two_norm2();
                                        return res;
                                    } ();

                    // add value to vector of interior neumann fluxes
                    facetCouplingFluxes[data.localIndexInInteractionVolume()] += getX(facetVolVars)*
                                                                                 curScvf.area()*
                                                                                 elemVolVars[curScvf.insideScvIdx()].extrusionFactor()*
                                                                                 MpfaHelper::nT_M_v(n, facetVolVars.diffusionCoefficient(phaseIdx, compIdx), v);
                }
                else
                {
                    // get the complete data of the actual interior neumann face
                    const auto completeFacetData = data.completeCoupledFacetData(fvGeometry);

                    // calculate "leakage factor"
                    const auto n = curScvf.unitOuterNormal();
                    const auto v = [&] ()
                                    {
                                        auto res = n;
                                        res *= -0.5*completeFacetData.volVars().extrusionFactor();
                                        res /= res.two_norm2();
                                        return res;
                                    } ();

                    auto getD = TensorLambdaFactory<TypeTag, myDiscretizationMethod>::getDiffusionLambda(phaseIdx, compIdx);
                    const auto D = getD(completeFacetData.problem(),
                                        completeFacetData.element(),
                                        completeFacetData.volVars(),
                                        completeFacetData.fvGeometry(),
                                        completeFacetData.scv());

                    // add value to vector of interior neumann fluxes
                    facetCouplingFluxes[data.localIndexInInteractionVolume()] += getX(completeFacetData.volVars())*
                                                                                 curScvf.area()*
                                                                                 elemVolVars[curScvf.insideScvIdx()].extrusionFactor()*
                                                                                 MpfaHelper::nT_M_v(n, D, v);
                }
            }
        }

        return flux + cij*facetCouplingFluxes;
    }
};

} // end namespace

#endif
