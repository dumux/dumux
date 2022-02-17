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
  * \ingroup PNMTwoPModel
  * \brief Flux variables cache for the two-phase-flow PNM
  */
#ifndef DUMUX_PNM_2P_FLUXVARIABLESCACHE_HH
#define DUMUX_PNM_2P_FLUXVARIABLESCACHE_HH

#include <array>
#include <algorithm>
#include <dune/common/reservedvector.hh>
#include <dumux/porenetwork/common/throatproperties.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMTwoPModel
 * \brief Flux variables cache for the two-phase-flow PNM
 *        Store data required for flux calculation.
 */
template<class AdvectionType, int maxNumCorners = 4>
class TwoPFluxVariablesCache
{
    using Scalar = typename AdvectionType::Scalar;
    static constexpr auto numPhases = 2;
    using NumCornerVector = Dune::ReservedVector<Scalar, maxNumCorners>;

public:
    TwoPFluxVariablesCache()
    {
        deltaPc_ = getParam<Scalar>("Regularization.DeltaPc");
        intervalPc_ = getParam<Scalar>("Regularization.IntervalPc");
    }

    template<class Problem, class Element, class FVElementGeometry,
             class ElementVolumeVariables, class SubControlVolumeFace>
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace& scvf,
                bool invaded)
    {
        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
        throatCrossSectionShape_ = fvGeometry.gridGeometry().throatCrossSectionShape(eIdx);
        throatShapeFactor_ = fvGeometry.gridGeometry().throatShapeFactor(eIdx);
        pc_ = std::max(elemVolVars[0].capillaryPressure(), elemVolVars[1].capillaryPressure());
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

        const auto& cornerHalfAngles = spatialParams.cornerHalfAngles(element);

        wettingLayerArea_.clear(); wettingLayerArea_.resize(cornerHalfAngles.size());
        const Scalar totalThroatCrossSectionalArea = spatialParams.throatCrossSectionalArea(element, elemVolVars);
        entryWettingLayerArea_.clear(); entryWettingLayerArea_.resize(cornerHalfAngles.size());
        snapoffWettingLayerArea_.clear(); snapoffWettingLayerArea_.resize(cornerHalfAngles.size());
        deltaPcEntryWettingLayerArea_.clear(); deltaPcEntryWettingLayerArea_.resize(cornerHalfAngles.size());
        deltaPcSnapoffWettingLayerArea_.clear(); deltaPcSnapoffWettingLayerArea_.resize(cornerHalfAngles.size());

        // interval pcEntry corresponding wetting layer
        intervalPcEntryWettingLayerArea_.clear(); intervalPcEntryWettingLayerArea_.resize(cornerHalfAngles.size());
        deltaIntervalPcEntryWettingLayerArea_.clear(), deltaIntervalPcEntryWettingLayerArea_.resize(cornerHalfAngles.size());

        const Scalar theta = spatialParams.contactAngle(element, elemVolVars);
        for (int i = 0; i< cornerHalfAngles.size(); ++i)
        {
            wettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(curvatureRadius(), theta, cornerHalfAngles[i]);
            entryWettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(entryCurvatureRadius(), theta, cornerHalfAngles[i]);
            deltaPcEntryWettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(delatEntryCurvatureRadius(), theta, cornerHalfAngles[i]);
            deltaPcSnapoffWettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(delatSnapoffCurvatureRadius(), theta, cornerHalfAngles[i]);
            snapoffWettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(snapoffCurvatureRadius(), theta, cornerHalfAngles[i]);
            intervalPcEntryWettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea((surfaceTension_/(pcEntry_ + intervalPc_)), theta, cornerHalfAngles[i]);
            deltaIntervalPcEntryWettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea((surfaceTension_/(pcEntry_ + intervalPc_ + deltaPc_)), theta, cornerHalfAngles[i]);
        }

        // make sure the wetting phase area does not exceed the total cross-section area
        throatCrossSectionalArea_[wPhaseIdx()] = std::min(
            std::accumulate(wettingLayerArea_.begin(), wettingLayerArea_.end(), 0.0),
            totalThroatCrossSectionalArea
        );
        throatCrossSectionalArea_[nPhaseIdx()] = totalThroatCrossSectionalArea - throatCrossSectionalArea_[wPhaseIdx()];

        // returns the throat cross sectional area for both phases at critical status corresponds entry pressure
        entryThroatCrossSectionalArea_[wPhaseIdx()] = std::min(std::accumulate(entryWettingLayerArea_.begin(), entryWettingLayerArea_.end(), 0.0), totalThroatCrossSectionalArea);
        entryThroatCrossSectionalArea_[nPhaseIdx()] = spatialParams.throatCrossSectionalArea(element, elemVolVars) - entryThroatCrossSectionalArea_[wPhaseIdx()];

        // returns the throat cross sectional area for both phases at critical status corresponds snapoff pressure
        snapoffThroatCrossSectionalArea_[wPhaseIdx()] = std::min(std::accumulate(snapoffWettingLayerArea_.begin(), snapoffWettingLayerArea_.end(), 0.0), totalThroatCrossSectionalArea);
        snapoffThroatCrossSectionalArea_[nPhaseIdx()] = spatialParams.throatCrossSectionalArea(element, elemVolVars) - snapoffThroatCrossSectionalArea_[wPhaseIdx()];

        auto deltaPcEntryThroatCrossSectionalArea = std::min(std::accumulate(deltaPcEntryWettingLayerArea_.begin(), deltaPcEntryWettingLayerArea_.end(), 0.0), totalThroatCrossSectionalArea);
        auto deltaPcSnapoffThroatCrossSectionalArea = std::min(std::accumulate(deltaPcSnapoffWettingLayerArea_.begin(), deltaPcSnapoffWettingLayerArea_.end(), 0.0), totalThroatCrossSectionalArea);
        auto intervalPcEntryThroatCrossSectionalArea =  std::min(std::accumulate(intervalPcEntryWettingLayerArea_.begin(), intervalPcEntryWettingLayerArea_.end(), 0.0), totalThroatCrossSectionalArea);
        auto deltaIntervalPcEntryThroatCrossSectionalArea = std::min(std::accumulate(deltaIntervalPcEntryWettingLayerArea_.begin(), deltaIntervalPcEntryWettingLayerArea_.end(), 0.0), totalThroatCrossSectionalArea);

        deltaNonWettingCrossSectionalAreaEntry_ = totalThroatCrossSectionalArea - deltaPcEntryThroatCrossSectionalArea;
        deltaNonWettingCrossSectionalAreaSnapoff_ = totalThroatCrossSectionalArea - deltaPcSnapoffThroatCrossSectionalArea;
        intervalNonWettingCrossSectionalEntry_ = totalThroatCrossSectionalArea - intervalPcEntryThroatCrossSectionalArea;
        deltaIntervalNonWettingCrossSectionalEntry_ = totalThroatCrossSectionalArea - deltaIntervalPcEntryThroatCrossSectionalArea;

        entryDerivativeCrossSectionalAreaPc_[wPhaseIdx()] = (deltaPcEntryThroatCrossSectionalArea - throatCrossSectionalArea_[wPhaseIdx()])/deltaPc_;
        entryDerivativeCrossSectionalAreaPc_[nPhaseIdx()] = - entryDerivativeCrossSectionalAreaPc_[wPhaseIdx()];
        snapoffDerivativeCrossSectionalAreaPc_[wPhaseIdx()] = (deltaPcSnapoffThroatCrossSectionalArea - throatCrossSectionalArea_[wPhaseIdx()])/deltaPc_;
        snapoffDerivativeCrossSectionalAreaPc_[wPhaseIdx()] = - snapoffDerivativeCrossSectionalAreaPc_[nPhaseIdx()];

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            singlePhaseCache_.fill(problem, element, fvGeometry, scvf, elemVolVars, *this, phaseIdx);
            nonWettingPhaseCache_.fill(problem, element, fvGeometry, scvf, elemVolVars, *this, phaseIdx);
            wettingLayerCache_.fill(problem, element, fvGeometry, scvf, elemVolVars, *this, phaseIdx);
        }

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            transmissibility_[phaseIdx] = AdvectionType::calculateTransmissibility(
                problem, element, fvGeometry, scvf, elemVolVars, *this, phaseIdx
            );
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
     * \brief Returns the curvature radius within the throat corresponds to entry pressure.
     */
    Scalar entryCurvatureRadius() const
    { return surfaceTension_ / pcEntry_;}

    Scalar delatEntryCurvatureRadius() const
    { return surfaceTension_ / (pcEntry_ + deltaPc_); }

    Scalar delatSnapoffCurvatureRadius() const
    { return surfaceTension_ / (pcSnapoff_ + deltaPc_); }

    /*!
     * \brief Returns the curvature radius within the throat corresponds to snap off pressure.
     */
    Scalar snapoffCurvatureRadius() const
    { return surfaceTension_ / pcSnapoff_;}

    /*!
     * \brief Returns the cross-sectional area of a wetting layer within
     *        one of the throat's corners.
     */
    Scalar wettingLayerCrossSectionalArea(const int cornerIdx) const
    { return wettingLayerArea_[cornerIdx]; }

    Scalar entryWettingLayerArea(const int cornerIdx) const
    { return entryWettingLayerArea_[cornerIdx]; }

    Scalar snapoffWettingLayerArea(const int cornerIdx) const
    { return snapoffWettingLayerArea_[cornerIdx]; }

    Scalar deltaPcSnapoffWettingLayerArea(const int cornerIdx) const
    { return deltaPcSnapoffWettingLayerArea_[cornerIdx]; }

    Scalar deltaPcEntryWettingLayerArea(const int cornerIdx) const
    { return deltaPcEntryWettingLayerArea_[cornerIdx]; }

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


    Scalar entryThroatCrossSectionalArea(const int phaseIdx) const
    { return entryThroatCrossSectionalArea_[phaseIdx]; }

    Scalar snapoffThroatCrossSectionalArea(const int phaseIdx) const
    { return snapoffThroatCrossSectionalArea_[phaseIdx]; }

    Scalar entryDerivativeCrossSectionalAreaPc(const int phaseIdx) const
    { return entryDerivativeCrossSectionalAreaPc_[phaseIdx]; }

    Scalar snapoffDerivativeCrossSectionalAreaPc(const int phaseIdx) const
    { return snapoffDerivativeCrossSectionalAreaPc_[phaseIdx]; }

    Scalar deltaPc() const
    { return deltaPc_; }

    Scalar deltaNonWettingCrossSectionalAreaEntry() const
    { return deltaNonWettingCrossSectionalAreaEntry_; }

    Scalar deltaNonWettingCrossSectionalAreaSnapoff() const
    { return deltaNonWettingCrossSectionalAreaSnapoff_; }

    Scalar intervalNonWettingCrossSectionalEntry() const
    { return intervalNonWettingCrossSectionalEntry_; }

    Scalar deltaIntervalNonWettingCrossSectionalEntry() const
    { return deltaIntervalNonWettingCrossSectionalEntry_; }

    Scalar intervalPc() const
    { return intervalPc_; }

private:
    Throat::Shape throatCrossSectionShape_;
    Scalar throatShapeFactor_;
    std::array<Scalar, numPhases> transmissibility_;
    std::array<Scalar, numPhases> throatCrossSectionalArea_;
    std::array<Scalar, numPhases> entryDerivativeCrossSectionalAreaPc_;
    std::array<Scalar, numPhases> snapoffDerivativeCrossSectionalAreaPc_;

    Scalar throatLength_;
    Scalar throatInscribedRadius_;
    Scalar pcEntry_;
    Scalar pcSnapoff_;
    Scalar pc_;
    Scalar deltaPc_;
    Scalar intervalPc_, intervalNonWettingCrossSectionalEntry_, deltaIntervalNonWettingCrossSectionalEntry_;
    Scalar surfaceTension_;
    bool invaded_;
    NumCornerVector wettingLayerArea_, entryWettingLayerArea_, snapoffWettingLayerArea_;
    NumCornerVector deltaPcEntryWettingLayerArea_, deltaPcSnapoffWettingLayerArea_;
    NumCornerVector intervalPcEntryWettingLayerArea_, deltaIntervalPcEntryWettingLayerArea_;

    std::size_t nPhaseIdx_;
    Scalar poreToPoreDistance_;

    typename AdvectionType::Transmissibility::SinglePhaseCache singlePhaseCache_;
    typename AdvectionType::Transmissibility::NonWettingPhaseCache nonWettingPhaseCache_;
    typename AdvectionType::Transmissibility::WettingLayerCache wettingLayerCache_;

    std::array<Scalar, numPhases> entryThroatCrossSectionalArea_;
    std::array<Scalar, numPhases> snapoffThroatCrossSectionalArea_;

    Scalar deltaNonWettingCrossSectionalAreaEntry_;
    Scalar deltaNonWettingCrossSectionalAreaSnapoff_;
};

} // end Dumux::PoreNetwork

#endif
