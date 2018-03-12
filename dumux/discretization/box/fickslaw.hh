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
 *        diffusive mass fluxes due to molecular diffusion with Fick's law.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_FICKS_LAW_HH

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{
// forward declaration
template<class TypeTag, DiscretizationMethod discMethod>
class FicksLawImplementation;

/*!
 * \ingroup BoxFicksLaw
 * \brief Specialization of Fick's Law for the box method.
 */
template <class TypeTag>
class FicksLawImplementation<TypeTag, DiscretizationMethod::box>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVarCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using BalanceEqOpts = typename GET_PROP_TYPE(TypeTag, BalanceEqOpts);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using Element = typename GridView::template Codim<0>::Entity;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum
    {
        numPhases = GET_PROP_TYPE(TypeTag, ModelTraits)::numPhases(),
        numComponents =  GET_PROP_TYPE(TypeTag, ModelTraits)::numComponents()
    };
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;

public:

    static ComponentFluxVector flux(const Problem& problem,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolumeFace& scvf,
                                    const int phaseIdx,
                                    const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        ComponentFluxVector componentFlux(0.0);

        // get inside and outside diffusion tensors and calculate the harmonic mean
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        // evaluate gradX at integration point and interpolate density
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto& shapeValues = fluxVarsCache.shapeValues();

        // density interpolation
        Scalar rho(0.0);
        for (auto&& scv : scvs(fvGeometry))
            rho += elemVolVars[scv].molarDensity(phaseIdx)*shapeValues[scv.indexInElement()][0];

        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            if(compIdx == FluidSystem::getMainComponent(phaseIdx))
                continue;

            // effective diffusion tensors
            using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);
            auto insideD = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(),
                                                              insideVolVars.saturation(phaseIdx),
                                                              insideVolVars.diffusionCoefficient(phaseIdx, compIdx));
            auto outsideD = EffDiffModel::effectiveDiffusivity(outsideVolVars.porosity(),
                                                               outsideVolVars.saturation(phaseIdx),
                                                               outsideVolVars.diffusionCoefficient(phaseIdx, compIdx));

            // scale by extrusion factor
            insideD *= insideVolVars.extrusionFactor();
            outsideD *= outsideVolVars.extrusionFactor();

            // the resulting averaged diffusion tensor
            const auto D = problem.spatialParams().harmonicMean(insideD, outsideD, scvf.unitOuterNormal());

            // the mole/mass fraction gradient
            GlobalPosition gradX(0.0);
            for (auto&& scv : scvs(fvGeometry))
                gradX.axpy(elemVolVars[scv].moleFraction(phaseIdx, compIdx), fluxVarsCache.gradN(scv.indexInElement()));

            // compute the diffusive flux
            componentFlux[compIdx] = -1.0*rho*vtmv(scvf.unitOuterNormal(), D, gradX)*scvf.area();
            if (BalanceEqOpts::mainComponentIsBalanced(phaseIdx) && !FluidSystem::isTracerFluidSystem())
                componentFlux[phaseIdx] -= componentFlux[compIdx];
        }
        return componentFlux;
    }

    // compute transmissibilities ti for analytical jacobians
    static std::array<std::vector<Scalar>, numComponents>
    calculateTransmissibilities(const Problem& problem,
                                const Element& element,
                                const FVElementGeometry& fvGeometry,
                                const ElementVolumeVariables& elemVolVars,
                                const SubControlVolumeFace& scvf,
                                const FluxVarCache& fluxVarCache,
                                const int phaseIdx)
    {
        Scalar rho(0.0);
        const auto& shapeValues = fluxVarCache.shapeValues();
        for (auto&& scv : scvs(fvGeometry))
            rho += elemVolVars[scv].molarDensity(phaseIdx)*shapeValues[scv.indexInElement()][0];

        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        std::array<std::vector<Scalar>, numComponents> ti;
        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            if(compIdx == FluidSystem::getMainComponent(phaseIdx))
                continue;

            // effective diffusion tensors
            using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);
            auto insideD = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(),
                                                              insideVolVars.saturation(phaseIdx),
                                                              insideVolVars.diffusionCoefficient(phaseIdx, compIdx));
            auto outsideD = EffDiffModel::effectiveDiffusivity(outsideVolVars.porosity(),
                                                               outsideVolVars.saturation(phaseIdx),
                                                               outsideVolVars.diffusionCoefficient(phaseIdx, compIdx));

            // scale by extrusion factor
            insideD *= insideVolVars.extrusionFactor();
            outsideD *= outsideVolVars.extrusionFactor();

            // the resulting averaged diffusion tensor
            const auto D = problem.spatialParams().harmonicMean(insideD, outsideD, scvf.unitOuterNormal());

            ti[compIdx].resize(fvGeometry.numScv());
            for (auto&& scv : scvs(fvGeometry))
                ti[compIdx][scv.indexInElement()] = -rho*vtmv(scvf.unitOuterNormal(), D, fluxVarCache.gradN(scv.indexInElement()))*scvf.area();
        }

        return ti;
    }
};

} // end namespace Dumux

#endif
