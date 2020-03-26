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
 * \ingroup BoxFlux
 * \brief This file contains the data which is required to calculate
 *        diffusive fluxes due to molecular diffusion with Fick's law.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_FICKS_LAW_HH

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/deprecated.hh>
#include <dumux/discretization/method.hh>

#include <dumux/flux/fickiandiffusioncoefficients.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, DiscretizationMethod discMethod, ReferenceSystemFormulation referenceSystem>
class FicksLawImplementation;

/*!
 * \ingroup BoxFlux
 * \brief Specialization of Fick's Law for the box method.
 */
template <class TypeTag, ReferenceSystemFormulation referenceSystem>
class FicksLawImplementation<TypeTag, DiscretizationMethod::box, referenceSystem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using FluxVarCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
    using BalanceEqOpts = GetPropType<TypeTag, Properties::BalanceEqOpts>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum
    {
        numPhases = ModelTraits::numFluidPhases(),
        numComponents = ModelTraits::numFluidComponents()
    };
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;

public:
    using DiffusionCoefficientsContainer = FickianDiffusionCoefficients<Scalar, numPhases, numComponents>;

    //return the reference system
    static constexpr ReferenceSystemFormulation referenceSystemFormulation()
    { return referenceSystem; }

    static ComponentFluxVector flux(const Problem& problem,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolumeFace& scvf,
                                    const int phaseIdx,
                                    const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        // get inside and outside diffusion tensors and calculate the harmonic mean
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        // evaluate gradX at integration point and interpolate density
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto& shapeValues = fluxVarsCache.shapeValues();

        // density interpolation
        Scalar rho(0.0);
        for (auto&& scv : scvs(fvGeometry))
            rho += massOrMolarDensity(elemVolVars[scv], referenceSystem, phaseIdx)*shapeValues[scv.indexInElement()][0];

        ComponentFluxVector componentFlux(0.0);
        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            if constexpr (!FluidSystem::isTracerFluidSystem())
                if (compIdx == FluidSystem::getMainComponent(phaseIdx))
                    continue;

            const auto D = averageDiffusionCoefficient_(phaseIdx, compIdx, insideVolVars, outsideVolVars, problem, scvf);

            // compute the diffusive flux
            const auto massOrMoleFrac = [&](const SubControlVolume& scv){ return massOrMoleFraction(elemVolVars[scv], referenceSystem, phaseIdx, compIdx); };
            componentFlux[compIdx] = discreteFlux_(fvGeometry, scvf, fluxVarsCache, massOrMoleFrac, D, rho);

            // if the main component is balanced subtract the same flux from there (conservation)
            if constexpr (!FluidSystem::isTracerFluidSystem())
                if (BalanceEqOpts::mainComponentIsBalanced(phaseIdx))
                    componentFlux[FluidSystem::getMainComponent(phaseIdx)] -= componentFlux[compIdx];
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
            rho += massOrMolarDensity(elemVolVars[scv], referenceSystem, phaseIdx)*shapeValues[scv.indexInElement()][0];

        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        std::array<std::vector<Scalar>, numComponents> ti;
        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            if constexpr (!FluidSystem::isTracerFluidSystem())
                if (compIdx == FluidSystem::getMainComponent(phaseIdx))
                    continue;

            const auto D = averageDiffusionCoefficient_(phaseIdx, compIdx, insideVolVars, outsideVolVars, problem, scvf);

            ti[compIdx].resize(fvGeometry.numScv());
            for (auto&& scv : scvs(fvGeometry))
                ti[compIdx][scv.indexInElement()] = -rho*vtmv(scvf.unitOuterNormal(), D, fluxVarCache.gradN(scv.indexInElement()))*scvf.area();
        }

        return ti;
    }

private:
    static Scalar averageDiffusionCoefficient_(const int phaseIdx, const int compIdx,
                                               const VolumeVariables& insideVV, const VolumeVariables& outsideVV,
                                               const Problem& problem,
                                               const SubControlVolumeFace& scvf)
    {
        // effective diffusion tensors
        auto [insideD, outsideD] = diffusionCoefficientsAtInterface_(phaseIdx, compIdx, insideVV, outsideVV);

        // scale by extrusion factor
        insideD *= insideVV.extrusionFactor();
        outsideD *= outsideVV.extrusionFactor();

        // the resulting averaged diffusion tensor
        return problem.spatialParams().harmonicMean(insideD, outsideD, scvf.unitOuterNormal());
    }

    static std::pair<Scalar, Scalar>
    diffusionCoefficientsAtInterface_(const int phaseIdx, const int compIdx,
                                      const VolumeVariables& insideVV, const VolumeVariables& outsideVV)
    {
        if constexpr (!FluidSystem::isTracerFluidSystem())
        {
            using EffDiffModel = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
            const auto mainCompIdx = FluidSystem::getMainComponent(phaseIdx);
            const auto insideD = Deprecated::template effectiveDiffusionCoefficient<EffDiffModel>(insideVV, phaseIdx, mainCompIdx, compIdx);
            const auto outsideD = Deprecated::template effectiveDiffusionCoefficient<EffDiffModel>(outsideVV, phaseIdx, mainCompIdx, compIdx);
            return { std::move(insideD), std::move(outsideD) };
        }
        else
        {
            using EffDiffModel = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
            const auto insideD = Deprecated::template effectiveDiffusionCoefficient<EffDiffModel>(insideVV, 0, 0, compIdx);
            const auto outsideD = Deprecated::template effectiveDiffusionCoefficient<EffDiffModel>(outsideVV, 0, 0, compIdx);
            return { std::move(insideD), std::move(outsideD) };
        }
    }

    template<class EvaluateVariable, class Tensor>
    static Scalar discreteFlux_(const FVElementGeometry& fvGeometry,
                                const SubControlVolumeFace& scvf,
                                const FluxVarCache& fluxVarsCache,
                                const EvaluateVariable& massOrMoleFraction,
                                const Tensor& D, const Scalar preFactor)
    {
        Dune::FieldVector<Scalar, dimWorld> gradX(0.0);
        for (auto&& scv : scvs(fvGeometry))
            gradX.axpy(massOrMoleFraction(scv), fluxVarsCache.gradN(scv.indexInElement()));
        return -1.0*preFactor*vtmv(scvf.unitOuterNormal(), D, gradX)*scvf.area();
    }
};

} // end namespace Dumux

#endif
