// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkFlux
 * \brief This file contains the data which is required to calculate
 *        the fluxes of the pore network model over a face of a finite volume.
 */
#ifndef DUMUX_FLUX_PNM_ADVECTION_HH
#define DUMUX_FLUX_PNM_ADVECTION_HH

#include <array>
#include <dumux/common/parameters.hh>

namespace Dumux::PoreNetwork::Detail {

template<class... TransmissibilityLawTypes>
struct Transmissibility : public TransmissibilityLawTypes... {};

} // end namespace Dumux::PoreNetwork::Detail

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PoreNetworkFlux
 * \brief Hagen–Poiseuille-type flux law to describe the advective flux for pore-network models.
 */
template<class ScalarT, class... TransmissibilityLawTypes>
class CreepingFlow
{

public:
    //! Export the Scalar type
    using Scalar = ScalarT;

    //! Export the transmissibility law
    using Transmissibility = Detail::Transmissibility<TransmissibilityLawTypes...>;

    /*!
     * \brief Returns the advective flux of a fluid phase
     *        across the given sub-control volume face (corresponding to a pore throat).
     * \note The flux is given in N*m, and can be converted
     *       into a volume flux (m^3/s) or mass flux (kg/s) by applying an upwind scheme
     *       for the mobility (1/viscosity) or the product of density and mobility, respectively.
     */
    template<class Problem, class Element, class FVElementGeometry,
             class ElementVolumeVariables, class SubControlVolumeFace, class ElemFluxVarsCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const int phaseIdx,
                       const ElemFluxVarsCache& elemFluxVarsCache)
    {
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        // Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];

        // calculate the pressure difference
        const Scalar deltaP = insideVolVars.pressure(phaseIdx) - outsideVolVars.pressure(phaseIdx);
        const Scalar transmissibility = fluxVarsCache.transmissibility(phaseIdx);
        using std::isfinite;    
        assert(isfinite(transmissibility));

        Scalar volumeFlow = transmissibility*deltaP;

        // add gravity term
        static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");
        if (enableGravity)
        {
            const Scalar rho = 0.5*insideVolVars.density(phaseIdx) + 0.5*outsideVolVars.density(phaseIdx);
            const Scalar g = problem.spatialParams().gravity(scvf.center()) * scvf.unitOuterNormal();

            // The transmissibility is with respect to the effective throat length (potentially dropping the pore body radii).
            // For gravity, we need to consider the total throat length (i.e., the cell-center to cell-center distance).
            // This might cause some inconsistencies TODO: is there a better way?
            volumeFlow += transmissibility * fluxVarsCache.poreToPoreDistance() * rho * g;
        }

        return volumeFlow;
    }

    /*!
     * \brief Returns the throat conductivity
     */
    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
    static Scalar calculateTransmissibility(const Problem& problem,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                            const ElementVolumeVariables& elemVolVars,
                                            const FluxVariablesCache& fluxVarsCache,
                                            const int phaseIdx)
    {
        if constexpr (ElementVolumeVariables::VolumeVariables::numFluidPhases() == 1)
            return Transmissibility::singlePhaseTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVarsCache, phaseIdx);
        else
        {
            static_assert(ElementVolumeVariables::VolumeVariables::numFluidPhases() == 2);

            const auto& spatialParams = problem.spatialParams();
            using FluidSystem = typename ElementVolumeVariables::VolumeVariables::FluidSystem;
            const int wPhaseIdx = spatialParams.template wettingPhase<FluidSystem>(element, elemVolVars);
            const bool invaded = fluxVarsCache.invaded();

            if (phaseIdx == wPhaseIdx)
            {
                return invaded ? Transmissibility::wettingLayerTransmissibility(element, fvGeometry, scvf, fluxVarsCache)
                               : Transmissibility::singlePhaseTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVarsCache, phaseIdx);
            }
            else // non-wetting phase
            {
                return invaded ? Transmissibility::nonWettingPhaseTransmissibility(element, fvGeometry, scvf, fluxVarsCache)
                               : 0.0;
            }
        }
    }

    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
    static std::array<Scalar, 2> calculateTransmissibilities(const Problem& problem,
                                                             const Element& element,
                                                             const FVElementGeometry& fvGeometry,
                                                             const ElementVolumeVariables& elemVolVars,
                                                             const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                                             const FluxVariablesCache& fluxVarsCache)
    {
        static_assert(ElementVolumeVariables::VolumeVariables::numFluidPhases() == 1);
        const Scalar t = calculateTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVarsCache, 0);
        return {t, -t};
    }
};

/*!
 * \ingroup PoreNetworkFlux
 * \brief Non-creeping flow flux law to describe the advective flux for pore-network models based on
 *        El-Zehairy et al.(2019).
 *        https://doi.org/10.1016/j.advwatres.2019.103378
 */
template <class ScalarT, class... TransmissibilityLawTypes>
class NonCreepingFlow
{
public:
    //! Export the Scalar type
    using Scalar = ScalarT;

    //! Export the creeping flow transmissibility law
    using TransmissibilityCreepingFlow = Detail::Transmissibility<TransmissibilityLawTypes...>;

    //! Inherit transmissibility from creeping flow transmissibility law to cache non-creeping flow-related parameters
    class Transmissibility: public TransmissibilityCreepingFlow
    {
    public:
        class SinglePhaseCache : public TransmissibilityCreepingFlow::SinglePhaseCache
        {
            using ParentType = typename TransmissibilityCreepingFlow::SinglePhaseCache;
        public:
            using Scalar = ScalarT;

            template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
            void fill(const Problem& problem,
                    const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const typename FVElementGeometry::SubControlVolumeFace& scvf,
                    const ElementVolumeVariables& elemVolVars,
                    const FluxVariablesCache& fluxVarsCache,
                    const int phaseIdx)
            {
                ParentType::fill(problem, element, fvGeometry,scvf, elemVolVars, fluxVarsCache, phaseIdx);
                const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);

                for (const auto& scv : scvs(fvGeometry))
                {
                    const auto localIdx = scv.indexInElement();
                    const Scalar throatToPoreAreaRatio = fluxVarsCache.throatCrossSectionalArea() / problem.spatialParams().poreCrossSectionalArea(element, scv, elemSol);

                    // dimensionless momentum coefficient
                    const Scalar momentumCoefficient = problem.spatialParams().momentumCoefficient(element, scv, elemSol);

                    // dimensionless kinetik-energy coefficient
                    const Scalar kineticEnergyCoefficient = problem.spatialParams().kineticEnergyCoefficient(element, scv, elemSol);

                    expansionCoefficient_[localIdx] = getExpansionCoefficient_(throatToPoreAreaRatio, momentumCoefficient, kineticEnergyCoefficient);
                    contractionCoefficient_[localIdx] = getContractionCoefficient_(throatToPoreAreaRatio, momentumCoefficient, kineticEnergyCoefficient);
                }
            }

            Scalar expansionCoefficient(int downstreamIdx) const
            { return expansionCoefficient_[downstreamIdx]; }

            Scalar contractionCoefficient(int upstreamIdx) const
            { return contractionCoefficient_[upstreamIdx]; }

        private:

            Scalar getExpansionCoefficient_(const Scalar throatToPoreAreaRatio, const Scalar momentumCoefficient, const Scalar kineticEnergyCoefficient) const
            {
                Scalar expansionCoefficient = throatToPoreAreaRatio * throatToPoreAreaRatio * (2 * momentumCoefficient - kineticEnergyCoefficient)
                                              + kineticEnergyCoefficient - 2 * momentumCoefficient * throatToPoreAreaRatio
                                              - (1 - throatToPoreAreaRatio * throatToPoreAreaRatio);

                return expansionCoefficient;
            }

            Scalar getContractionCoefficient_(const Scalar throatToPoreAreaRatio, const Scalar momentumCoefficient, const Scalar kineticEnergyCoefficient) const
            {
                const Scalar contractionAreaRatio = getContractionAreaRatio_(throatToPoreAreaRatio);
                Scalar contractionCoefficient = (1 - (throatToPoreAreaRatio * throatToPoreAreaRatio * kineticEnergyCoefficient - 2 * momentumCoefficient /*+1-throatToPoreAreaRatio*throatToPoreAreaRatio*/)
                                                     * contractionAreaRatio * contractionAreaRatio - 2 * contractionAreaRatio) / (contractionAreaRatio * contractionAreaRatio);

                return contractionCoefficient;
            }

            Scalar getContractionAreaRatio_(const Scalar throatToPoreAreaRatio) const
            {
                return 1.0-(1.0-throatToPoreAreaRatio)/(2.08*(1.0-throatToPoreAreaRatio)+0.5371);
            }

            std::array<Scalar, 2> expansionCoefficient_;
            std::array<Scalar, 2> contractionCoefficient_;
        };
    };

    /*!
     * \brief Returns the advective flux of a fluid phase
     *        across the given sub-control volume face (corresponding to a pore throat).
     * \note The flux is given in N*m, and can be converted
     *       into a volume flux (m^3/s) or mass flux (kg/s) by applying an upwind scheme
     *       for the mobility (1/viscosity) or the product of density and mobility, respectively.
     */
    template<class Problem, class Element, class FVElementGeometry,
             class ElementVolumeVariables, class SubControlVolumeFace, class ElemFluxVarsCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const int phaseIdx,
                       const ElemFluxVarsCache& elemFluxVarsCache)
    {
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        // Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];

        // calculate the pressure difference
        Scalar deltaP = insideVolVars.pressure(phaseIdx) - outsideVolVars.pressure(phaseIdx);

        // add gravity term
        static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");
        if (enableGravity)
        {
            const Scalar rho = 0.5*insideVolVars.density(phaseIdx) + 0.5*outsideVolVars.density(phaseIdx);

            // projected distance between pores in gravity field direction
            const auto poreToPoreVector = fluxVarsCache.poreToPoreDistance() * scvf.unitOuterNormal();
            const auto& gravityVector = problem.spatialParams().gravity(scvf.center());

            deltaP += rho * (gravityVector * poreToPoreVector);
        }
        const Scalar creepingFlowTransmissibility = fluxVarsCache.transmissibility(phaseIdx);
        const Scalar throatCrossSectionalArea = fluxVarsCache.throatCrossSectionalArea();
        using std::isfinite;
        assert(isfinite(creepingFlowTransmissibility));

        assert(scvf.insideScvIdx() == 0);
        assert(scvf.outsideScvIdx() == 1);

        // determine the flow direction to predict contraction and expansion
        const auto [upstreamIdx, downstreamIdx] = deltaP > 0 ? std::pair(0, 1) : std::pair(1, 0);

        const Scalar contractionCoefficient = fluxVarsCache.singlePhaseFlowVariables().contractionCoefficient(upstreamIdx);
        const Scalar expansionCoefficient = fluxVarsCache.singlePhaseFlowVariables().expansionCoefficient(downstreamIdx);
        const Scalar mu = elemVolVars[upstreamIdx].viscosity();

        //! El-Zehairy et al.(2019): Eq.(14):
        //! A0Q^2 + B0Q + C0 = 0  where Q is volumetric flowrate, A0 = [Kc + Ke] * rho /[2aij^2], Kc and Ke are contraction and expansion coefficient respectively,
        //! aij is cross sectional area of the throat, B0 = 1/K0 (K0 is trnasmissibility used in paper) and C0 = -deltaP
        //! In our implementation viscosity of fluid will be included using upwinding later. Thus, creepingFlowTransmissibility = mu * K0 and
        //! the volumetric flowrate calculated here q = mu * Q. Substitution of new variables into El-Zehairy et al.(2019), Eq.(14) gives:
        //! Aq^2 + Bq + C = 0  where
        //! A = A0 / mu^2,  B = B0/mu = 1/creepingFlowTransmissibility and C = C0
        //! attention: the q, volumetric flowrate, calculated here is always positive and its sign needs to be determined based on flow direction
        //! this approach is taken to prevent the term under the square root (discriminant) becoming negative
        const Scalar A = (contractionCoefficient * elemVolVars[upstreamIdx].density() + expansionCoefficient * elemVolVars[downstreamIdx].density())
                          / (2.0 * mu * mu * throatCrossSectionalArea * throatCrossSectionalArea);
        const Scalar B =  1.0/creepingFlowTransmissibility;
        using std::abs;
        const Scalar C = -abs(deltaP);

        using std::sqrt;
        const auto discriminant = B*B - 4*A*C;
        const auto q = (-B + sqrt(discriminant)) / (2*A);

        //! give the volume flowrate proper sign based on flow direction.
        if (upstreamIdx == 0)
            return q;
        else
            return -q;

    }

    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
    static Scalar calculateTransmissibility(const Problem& problem,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                            const ElementVolumeVariables& elemVolVars,
                                            const FluxVariablesCache& fluxVarsCache,
                                            const int phaseIdx)
    {
        static_assert(ElementVolumeVariables::VolumeVariables::numFluidPhases() == 1);
        return Transmissibility::singlePhaseTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVarsCache, phaseIdx);
    }

    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
    static std::array<Scalar, 2> calculateTransmissibilities(const Problem& problem,
                                                             const Element& element,
                                                             const FVElementGeometry& fvGeometry,
                                                             const ElementVolumeVariables& elemVolVars,
                                                             const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                                             const FluxVariablesCache& fluxVarsCache)
    {
        static_assert(ElementVolumeVariables::VolumeVariables::numFluidPhases() == 1);
        const Scalar t = calculateTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, fluxVarsCache, 0);
        return {t, -t};
    }

};

} // end namespace Dumux

#endif
