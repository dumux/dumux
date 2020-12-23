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
 * \brief This file contains the data which is required to calculate
 *        the fluxes of the pore network model over a face of a finite volume.
 */
#ifndef DUMUX_FLUX_PNM_ADVECTION_HH
#define DUMUX_FLUX_PNM_ADVECTION_HH

namespace Dumux
{

namespace Detail {

template<class... TransmissibilityLawTypes>
struct Transmissibility : public TransmissibilityLawTypes...
{};

}

/*!
 * \file
 * \ingroup PoreNetworkFlux
 * \brief Hagen–Poiseuille-type flux law to describe the advective flux for pore-network models.
 */
template<class ScalarT, class... TransmissibilityLawTypes>
class PoreNetworkCreepingFlow
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
        return std::array<Scalar, 2>{t, -t};
    }
};

/*!
 * \file
 * \ingroup PoreNetworkFlux
 * \brief Non-creeping flow flux law to describe the advective flux for pore-network models based on El-Zehairy et al.(2019).
 */
template <class ScalarT, class... TransmissibilityLawTypes>
class PoreNetworkNonCreepingFlow
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
                return 1-(1-throatToPoreAreaRatio)/(2.08*(1-throatToPoreAreaRatio)+0.5371);
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
        const Scalar deltaP = insideVolVars.pressure(phaseIdx) - outsideVolVars.pressure(phaseIdx);
        const Scalar creepingFlowTransmissibility = fluxVarsCache.transmissibility(phaseIdx);
        const Scalar throatCrossSectionalArea = fluxVarsCache.throatCrossSectionalArea();

        assert(scvf.insideScvIdx() == 0);
        assert(scvf.outsideScvIdx() == 1);

        // determine the flow direction to predict contraction and expansion
        const auto [upstreamIdx, downstreamIdx] = deltaP > 0 ? std::pair(0, 1) : std::pair(1, 0);

        const Scalar contractionCoefficient = fluxVarsCache.singlePhaseFlowVariables().contractionCoefficient(upstreamIdx);
        const Scalar expansionCoefficient = fluxVarsCache.singlePhaseFlowVariables().expansionCoefficient(downstreamIdx);
        const Scalar mu = elemVolVars[upstreamIdx].viscosity();

        //! A0q^2 + B0q + C0 = 0
        //! attention: the q, volumetric flowrate, calculated here is always positive and its sign need to be determined based on flow direction
        //! this approach is taken to prevent the term under the square root becoming negative
        const Scalar A0 = (contractionCoefficient * elemVolVars[upstreamIdx].density() + expansionCoefficient * elemVolVars[downstreamIdx].density())
                          / (2.0 * throatCrossSectionalArea * throatCrossSectionalArea);
        const Scalar B0 = mu / creepingFlowTransmissibility;
        const Scalar C0 = (upstreamIdx == 0) ? -deltaP: deltaP;

        using std::sqrt;
        const auto tmp0 = B0*B0 - 4*A0*C0;
        const auto q = (-B0 + sqrt(tmp0)) / (2*A0);

        //! give the volume flowrate proper sign based on flow direction
        if (upstreamIdx == 0)
            return mu * q;
        else
            return -mu * q;

        //TODO: gravity

    }

    /*!
     * \brief Returns the throat conductivity
     *
     * \param problem The problem
     * \param element The element
     * \param ElementVolumeVariables The element volume variables
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
        return std::array<Scalar, 2>{t, -t};
    }

};

} // end namespace

#endif
