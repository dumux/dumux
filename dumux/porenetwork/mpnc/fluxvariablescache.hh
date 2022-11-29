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
  * \ingroup PNMMPNCModel
  * \brief Flux variables cache for the MPNC PNM
  */
#ifndef DUMUX_PNM_MPNC_FLUXVARIABLESCACHE_HH
#define DUMUX_PNM_MPNC_FLUXVARIABLESCACHE_HH

#include <array>
#include <algorithm>
#include <dune/common/reservedvector.hh>
#include <dumux/porenetwork/common/throatproperties.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMMPNCModel
 * \brief Flux variables cache for the MPNC PNM
 *        Store data required for flux calculation.
 */
template<class AdvectionType, int maxNumCorners = 4>
class MPNCFluxVariablesCache
{
    using Scalar = typename AdvectionType::Scalar;
    static constexpr auto numPhases = 2;
    using NumCornerVector = Dune::ReservedVector<Scalar, maxNumCorners>;

public:

    template<class Problem, class Element, class FVElementGeometry,
             class ElementVolumeVariables, class SubControlVolumeFace>
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace& scvf,
                bool invaded,
                const Scalar factor)
    {
        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
        throatCrossSectionShape_ = fvGeometry.gridGeometry().throatCrossSectionShape(eIdx);
        throatShapeFactor_ = fvGeometry.gridGeometry().throatShapeFactor(eIdx);
        pcEntry_ = problem.spatialParams().pcEntry(element, elemVolVars);
        pcSnapoff_ = problem.spatialParams().pcSnapoff(element, elemVolVars);
        throatInscribedRadius_ = problem.spatialParams().throatInscribedRadius(element, elemVolVars);
        throatLength_ = problem.spatialParams().throatLength(element, elemVolVars);
        invaded_ = invaded;
        poreToPoreDistance_ = element.geometry().volume();

        // get the non-wetting phase index
        using FluidSystem = typename ElementVolumeVariables::VolumeVariables::FluidSystem;
        const auto& spatialParams = problem.spatialParams();
        nPhaseIdx_ = 1 - spatialParams.template wettingPhase<FluidSystem>(element, elemVolVars);

        // take the average surface tension of both adjacent pores TODO: is this correct?
        surfaceTension_ = 0.5*(elemVolVars[0].surfaceTension() + elemVolVars[1].surfaceTension());

        pcMax_ = pc_ = std::max(elemVolVars[0].capillaryPressure(), elemVolVars[1].capillaryPressure());
        pcAvg_ = 0.5 * (pc_ + std::min(elemVolVars[0].capillaryPressure(), elemVolVars[1].capillaryPressure()));

        if (!invaded)
            pc_ = std::max(elemVolVars[0].capillaryPressure(), elemVolVars[1].capillaryPressure());
        else
            pc_ = std::min(elemVolVars[0].capillaryPressure(), elemVolVars[1].capillaryPressure());

#if REGULARIZATIONWITHPRESSURE
        regInvasionInterval_[1] = pcEntry_ * (1 + regPercent_);
        regSnapoffInterval_[0] = pcSnapoff_* (1 - regPercent_);
#else
        const auto swEntry = fluidMatrixInteraction.sw(pcEntry_);
        const auto swRegEntry = swEntry - regPercent_;
        regInvasionInterval_[1] = fluidMatrixInteraction.pc(swRegEntry);

        const auto swSnapoff = fluidMatrixInteraction.sw(pcSnapoff_);
        const auto swRegSnapoff = swSnapoff + regPercent_;
        regSnapoffInterval_[0] = fluidMatrixInteraction.pc(swRegSnapoff);
#endif
        regInvasionInterval_[2] = regInvasionInterval_[1] + epsilonPc_;
        regSnapoffInterval_[2] = pcSnapoff_ + epsilonPc_;

        curvatureRadiusInvasion_[0] = surfaceTension_/regInvasionInterval_[1];
        curvatureRadiusInvasion_[1] = surfaceTension_/regInvasionInterval_[2];

        curvatureRadiusSnapoff_[0] = surfaceTension_/regSnapoffInterval_[1];
        curvatureRadiusSnapoff_[1] = surfaceTension_/regSnapoffInterval_[2];



        const auto& cornerHalfAngles = spatialParams.cornerHalfAngles(element);
        wettingLayerArea_.clear(); wettingLayerArea_.resize(cornerHalfAngles.size());

        if (invaded) // two-phase flow
        {
            const Scalar theta = spatialParams.contactAngle(element, elemVolVars);
            for (int i = 0; i< cornerHalfAngles.size(); ++i)
                wettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(curvatureRadius(), theta, cornerHalfAngles[i]);

            throatCrossSectionalArea_[wPhaseIdx()] = std::accumulate(wettingLayerArea_.begin(), wettingLayerArea_.end(), 0.0);
            throatCrossSectionalArea_[nPhaseIdx()] = spatialParams.throatCrossSectionalArea(element, elemVolVars) - throatCrossSectionalArea_[wPhaseIdx()];
        }
        else // single-phase flow
        {
            for (int i = 0; i< cornerHalfAngles.size(); ++i)
                wettingLayerArea_[i] = 0.0;

            throatCrossSectionalArea_[wPhaseIdx()] = spatialParams.throatCrossSectionalArea(element, elemVolVars);
            throatCrossSectionalArea_[nPhaseIdx()] = 0.0;
        }

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            singlePhaseCache_.fill(problem, element, fvGeometry, scvf, elemVolVars, *this, phaseIdx);
            nonWettingPhaseCache_.fill(problem, element, fvGeometry, scvf, elemVolVars, *this, phaseIdx);
            wettingLayerCache_.fill(problem, element, fvGeometry, scvf, elemVolVars, *this, phaseIdx);
        }

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            transmissibility_[phaseIdx] = AdvectionType::calculateTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, *this, phaseIdx) * factor;

        const auto filmThickness = [&](const Scalar pc, const Scalar beta)
        {
            using std::cos; using std::sin; using std::sqrt;
            const Scalar theta = spatialParams.contactAngle(element, elemVolVars);
            const Scalar rAm = this->surfaceTension() / pc;
            const Scalar b = rAm * cos(theta + beta) / sin(beta);
            return sqrt(b*b + rAm*rAm) - rAm;
        };

        // allow the breakdown of the film at some point
        if (invaded && cornerHalfAngles.size() > 0)
        {
            static const Scalar thresholdFactor = getParamFromGroup<Scalar>(problem.paramGroup(), "Problem.ThresholdFilmBreakdownFactor", 0.0);
            if (factor <= 0.0)
                return;

            bool breakUpOccurred = false;

            for (int i = 0; i< cornerHalfAngles.size(); ++i)
            {
                const Scalar wMax = filmThickness(this->pcEntry(), cornerHalfAngles[i]);
                const Scalar w = filmThickness(this->pc(), cornerHalfAngles[i]);

                if (w < thresholdFactor * wMax)
                {
                    std::cout << "Corner flow film breaks at " << eIdx << " (corner " << i << "). Max film thickness " << wMax << ", current thickness " << w << std::endl;
                    wettingLayerArea_[i] = 0.0;
                    breakUpOccurred = true;
                }
            }

            if (!breakUpOccurred)
                return;

            std::cout << "wPhase area before: " << throatCrossSectionalArea_[wPhaseIdx()] << ", nPhase area before: " << throatCrossSectionalArea_[nPhaseIdx()] << std::endl;
            throatCrossSectionalArea_[wPhaseIdx()] = std::accumulate(wettingLayerArea_.begin(), wettingLayerArea_.end(), 0.0);
            throatCrossSectionalArea_[nPhaseIdx()] = spatialParams.throatCrossSectionalArea(element, elemVolVars) - throatCrossSectionalArea_[wPhaseIdx()];
            std::cout << "wPhase area after: " << throatCrossSectionalArea_[wPhaseIdx()] << ", nPhase area after: " << throatCrossSectionalArea_[nPhaseIdx()] << std::endl;

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                singlePhaseCache_.fill(problem, element, fvGeometry, scvf, elemVolVars, *this, phaseIdx);
                nonWettingPhaseCache_.fill(problem, element, fvGeometry, scvf, elemVolVars, *this, phaseIdx);
                wettingLayerCache_.fill(problem, element, fvGeometry, scvf, elemVolVars, *this, phaseIdx);
            }

            std::cout << "wPhase transmissibility before: " << transmissibility_[wPhaseIdx()] << ", nPhase transmissibility before: " << transmissibility_[nPhaseIdx()] << std::endl;

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                transmissibility_[phaseIdx] = AdvectionType::calculateTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, *this, phaseIdx) * factor;

            std::cout << "wPhase transmissibility after: " << transmissibility_[wPhaseIdx()] << ", nPhase transmissibility after: " << transmissibility_[nPhaseIdx()] << std::endl;
        }
    }

    /*!
     * \brief Returns the throats's cross-sectional shape.
     */
    Throat::Shape throatCrossSectionShape() const
    { return throatCrossSectionShape_; }

    /*!
     * \brief Returns the throats's shape factor.
     */
    Scalar throatShapeFactor() const
    { return throatShapeFactor_; }

    /*!
     * \brief Returns the throats's transmissibility.
     */
    Scalar transmissibility(const int phaseIdx) const
    { return transmissibility_[phaseIdx]; }

    /*!
     * \brief Returns the throats's cross-sectional area for a given phaseIdx.
     */
    Scalar throatCrossSectionalArea(const int phaseIdx) const
    { return throatCrossSectionalArea_[phaseIdx]; }

    /*!
     * \brief Returns the throats's total cross-sectional area.
     */
    Scalar throatCrossSectionalArea() const
    { return throatCrossSectionalArea_[0] + throatCrossSectionalArea_[1]; }

    /*!
     * \brief Returns the throats's length.
     */
    Scalar throatLength() const
    { return throatLength_; }

    /*!
     * \brief Returns the throats's inscribed radius.
     */
    Scalar throatInscribedRadius() const
    { return throatInscribedRadius_; }

    /*!
     * \brief Returns the throats's entry capillary pressure.
     */
    Scalar pcEntry() const
    { return pcEntry_; }

    /*!
     * \brief Returns the throats's snap-off capillary pressure.
     */
    Scalar pcSnapoff() const
    { return pcSnapoff_; }

    /*!
     * \brief Returns the capillary pressure within the throat.
     */
    Scalar pc() const
    { return pc_; }

    /*!
     * \brief Returns the surface tension within the throat.
     */
    Scalar surfaceTension() const
    { return surfaceTension_; }

    /*!
     * \brief Returns true if the throat is invaded by the nonwetting phase.
     */
    bool invaded() const
    { return invaded_; }

    /*!
     * \brief Returns the curvature radius within the throat.
     */
    Scalar curvatureRadius() const
    { return surfaceTension_ / pc_;}

    /*!
     * \brief Returns the cross-sectional area of a wetting layer within
     *        one of the throat's corners.
     */
    Scalar wettingLayerCrossSectionalArea(const int cornerIdx) const
    { return wettingLayerArea_[cornerIdx]; }

    /*!
     * \brief Returns the index of the wetting phase.
     */
    std::size_t wPhaseIdx() const
    { return 1 - nPhaseIdx_; }

    /*!
     * \brief Returns the index of the nonwetting phase.
     */
    std::size_t nPhaseIdx() const
    { return nPhaseIdx_; }

    /*!
     * \brief Returns the throats's cached flow variables for single-phase flow.
     */
    const auto& singlePhaseFlowVariables() const
    { return singlePhaseCache_; }

    /*!
     * \brief Returns the throats's cached flow variables for the nonwetting phase.
     */
    const auto& nonWettingPhaseFlowVariables() const
    { return nonWettingPhaseCache_; }

    /*!
     * \brief Returns the throats's cached flow variables for the wetting phase.
     */
    const auto& wettingLayerFlowVariables() const
    { return wettingLayerCache_; }

    /*!
     * \brief Returns the throats's pore-to-pore-center distance.
     */
    Scalar poreToPoreDistance() const
    { return poreToPoreDistance_; }

private:
    Throat::Shape throatCrossSectionShape_;
    Scalar throatShapeFactor_;
    std::array<Scalar, numPhases> transmissibility_;
    std::array<Scalar, numPhases> throatCrossSectionalArea_;
    Scalar throatLength_;
    Scalar throatInscribedRadius_;
    Scalar pcEntry_;
    Scalar pcSnapoff_;
    Scalar pc_;
    Scalar surfaceTension_;
    bool invaded_;
    NumCornerVector wettingLayerArea_;
    std::size_t nPhaseIdx_;
    Scalar poreToPoreDistance_;

    typename AdvectionType::Transmissibility::SinglePhaseCache singlePhaseCache_;
    typename AdvectionType::Transmissibility::NonWettingPhaseCache nonWettingPhaseCache_;
    typename AdvectionType::Transmissibility::WettingLayerCache wettingLayerCache_;
};

} // end Dumux::PoreNetwork

#endif
