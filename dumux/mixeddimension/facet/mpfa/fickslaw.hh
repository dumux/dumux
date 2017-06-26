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
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

    // Always use the dynamic type for vectors (compatibility with the boundary)
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using CoefficientVector = typename BoundaryInteractionVolume::Traits::Vector;

public:

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
            // add neumann contributions
            if (data.faceType() == MpfaFaceTypes::interiorNeumann)
            {
                // get the scvf corresponding to actual interior neumann face
                const auto& curScvf = fvGeometry.scvf(data.scvfIndex());

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

                // get the diffusion coefficient from the lamda factory
                using LambdaFactory = TensorLambdaFactory<TypeTag, DiscretizationMethods::CCMpfa>;
                const auto D = LambdaFactory::getDiffusionLambda(phaseIdx, compIdx)(completeFacetData.problem(),
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
            // Add additional Dirichlet fluxes for interior Dirichlet faces
            else if (data.faceType() == MpfaFaceTypes::interiorDirichlet)
                flux += tij[startIdx + data.localIndexInInteractionVolume()]*getX(data.facetVolVars(fvGeometry));
        }

        return flux + cij*facetCouplingFluxes;
    }
};

} // end namespace

#endif
