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
 *        of Fick's Law. Specializations are provided for the different discretization methods.
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
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

public:
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

        // get the scaling factor for the effective diffusive fluxes
        auto factor = calculateEffectiveDiffusivityFactor(elemVolVars, scvf, phaseIdx);

        // if factor is zero, the flux will end up zero anyway
        if (factor == 0.0)
            return 0.0;

        // lambda functions depending on if we use mole or mass fractions
        auto xDensity = [useMoles, phaseIdx] (const VolumeVariables& volVars)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        auto xFraction = [useMoles, phaseIdx, compIdx] (const VolumeVariables& volVars)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        // calculate Tij*xj
        Scalar flux(0.0);
        unsigned int localIdx = 0;
        for (const auto volVarIdx : volVarsStencil)
            flux += tij[localIdx++]*xFraction(elemVolVars[volVarIdx]);

        // return effective mass flux
        return flux*interpolateDensity(elemVolVars, scvf, xDensity)*factor;
    }

    static Stencil stencil(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const SubControlVolumeFace& scvf)
    {
        const auto& globalFvGeometry = problem.model().globalFvGeometry();

        // return the scv (element) indices in the interaction region
        if (globalFvGeometry.scvfTouchesBoundary(scvf))
            return globalFvGeometry.boundaryInteractionVolumeSeed(scvf).globalScvIndices();
        else
            return globalFvGeometry.interactionVolumeSeed(scvf).globalScvIndices();
    }

private:
    template<typename DensityFunction>
    static Scalar interpolateDensity(const ElementVolumeVariables& elemVolVars,
                                     const SubControlVolumeFace& scvf,
                                     const DensityFunction& rhoFunction)
    {
        // use arithmetic mean of the densities around the scvf
        if (!scvf.boundary())
        {
            Scalar rho = rhoFunction(elemVolVars[scvf.insideScvIdx()]);
            for (auto outsideIdx : scvf.outsideScvIndices())
                rho += rhoFunction(elemVolVars[outsideIdx]);
            return rho/(scvf.outsideScvIndices().size()+1);
        }
        else
            return rhoFunction(elemVolVars[scvf.outsideScvIdx()]);
    }

    //! Here we want to calculate the factors with which the diffusion coefficient has to be
    //! scaled to get the effective diffusivity. For this we use the effective diffusivity with
    //! a diffusion coefficient of 1.0 as input. Then we scale the transmissibilites during flux
    //! calculation (above) with the harmonic average of the two factors
    static Scalar calculateEffectiveDiffusivityFactor(const ElementVolumeVariables& elemVolVars,
                                                      const SubControlVolumeFace& scvf,
                                                      const unsigned int phaseIdx)
    {
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        auto factor = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(),
                                                         insideVolVars.saturation(phaseIdx),
                                                         /*Diffusion coefficient*/ 1.0);

        if (!scvf.boundary())
        {
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
            auto outsideFactor = EffDiffModel::effectiveDiffusivity(outsideVolVars.porosity(),
                                                                    outsideVolVars.saturation(phaseIdx),
                                                                    /*Diffusion coefficient*/ 1.0);

            // use the harmonic mean of the two
            factor = harmonicMean(factor, outsideFactor);
        }

        return factor;
    }
};

} // end namespace

#endif
