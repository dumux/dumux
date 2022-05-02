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

template<class AdvectionType, int maxNumCorners = 4, bool regularized = true>
class TwoPFluxVariablesCache;
/*!
 * \ingroup PNMTwoPModel
 * \brief Flux variables cache for the two-phase-flow PNM
 *        Store data required for flux calculation.
 */
template<class AdvectionType>
class TwoPFluxVariablesCache<AdvectionType, 4, true>
{
    using Scalar = typename AdvectionType::Scalar;
    static constexpr auto numPhases = 2;
    using NumCornerVector = Dune::ReservedVector<Scalar, 4>;

public:
    TwoPFluxVariablesCache()
    {
        deltaPc_ = getParam<Scalar>("Regularization.DeltaPc", 1e-8);
        regularizeWithSaturation_ = getParam<bool>("Regularization.RegularizeWithSaturation", false);
        regularPcInterval_ = getParam<Scalar>("Regularization.IntervalPc", 1e2);
        regularizationPosition_ = getParam<Scalar>("Regularization.Position", 1.0);
        regSaturationPercentage_ = getParam<Scalar>("Regularization.SwPercentage", 1.0);
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

        // if (!regularizeWithSaturation_)
        // {
            // calculate boundary pressure values for invasion event
            regularizationPcEntry_[0] = pcEntry_ + regularPcInterval_*(regularizationPosition_-1);
            regularizationPcEntry_[1] = pcEntry_ + regularPcInterval_*regularizationPosition_;

            // calculate boundary pressure values for snap-off event
            regularizationPcSnapoff_[0] = pcSnapoff_ - regularPcInterval_*regularizationPosition_;
            regularizationPcSnapoff_[1] = pcSnapoff_ + regularPcInterval_*(1 - regularizationPosition_);
        // }
        // else
        // {
        //     const auto& spatialParams = problem.spatialParams();
        //     const auto fluidMatrixInteraction = spatialParams.fluidMatrixInteraction(element, scv, elemSol);
        //     const auto swEntry = fluidMatrixInteraction.sw(pcEntry_);
        //     const auto swRegEntry = swEntry - regSaturationPercentage_;
        //     regularizationPcEntry_[0] = pcEntry_;
        //     regularizationPcEntry_[1] = fluidMatrixInteraction.pc(swRegEntry);
        //     const auto swSnapoff = fluidMatrixInteraction.sw(pcSnapoff_);
        //     const auto swRegSnapoff = swSnapoff + regSaturationPercentage_;

        //     regularizationPcSnapoff_[0] = fluidMatrixInteraction.pc(swRegSnapoff);
        //     regularizationPcSnapoff_[1] = pcSnapoff_;
        // }


        regularizationPcEntry_[2] = regularizationPcEntry_[1] + deltaPc_;
        regularizationPcSnapoff_[2] = regularizationPcSnapoff_[1] + deltaPc_;

        // curvature at right boundary
        curvatureRadiusEntry_[0] = surfaceTension_/regularizationPcEntry_[1];
        curvatureRadiusEntry_[1] = surfaceTension_/regularizationPcEntry_[2];
        curvatureRadiusSnapoff_[0] = surfaceTension_/regularizationPcSnapoff_[1];
        curvatureRadiusSnapoff_[1] = surfaceTension_/regularizationPcSnapoff_[2];

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

        const Scalar theta = spatialParams.contactAngle(element, elemVolVars);
        for (int i = 0; i< cornerHalfAngles.size(); ++i)
        {
            wettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(curvatureRadius(), theta, cornerHalfAngles[i]);
            entryWettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(curvatureRadiusEntry(0), theta, cornerHalfAngles[i]);
            deltaPcEntryWettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(curvatureRadiusEntry(1), theta, cornerHalfAngles[i]);
            snapoffWettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(curvatureRadiusSnapoff(0), theta, cornerHalfAngles[i]);
            deltaPcSnapoffWettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(curvatureRadiusSnapoff(1), theta, cornerHalfAngles[i]);
        }

        // make sure the wetting phase area does not exceed the total cross-section area
        throatCrossSectionalArea_[wPhaseIdx()] = std::min(
            std::accumulate(wettingLayerArea_.begin(), wettingLayerArea_.end(), 0.0),
            totalThroatCrossSectionalArea
        );
        throatCrossSectionalArea_[nPhaseIdx()] = totalThroatCrossSectionalArea - throatCrossSectionalArea_[wPhaseIdx()];

        // return wetting area for regularization boundary for invasion
        regularBoundaryWettingThroatAreaEntry_[0] = std::min(std::accumulate(entryWettingLayerArea_.begin(), entryWettingLayerArea_.end(), 0.0), totalThroatCrossSectionalArea);
        regularBoundaryWettingThroatAreaEntry_[1] = std::min(std::accumulate(deltaPcEntryWettingLayerArea_.begin(), deltaPcEntryWettingLayerArea_.end(), 0.0), totalThroatCrossSectionalArea);
        regularBoundaryNonWettingThroatAreaEntry_[0] = totalThroatCrossSectionalArea - regularBoundaryWettingThroatAreaEntry_[0];
        regularBoundaryNonWettingThroatAreaEntry_[1] = totalThroatCrossSectionalArea - regularBoundaryWettingThroatAreaEntry_[1];

        // returns non wetting throat area for regularization boundary for snap-off
        regularBoundaryWettingThroatAreaSnapoff_[0] = std::min(std::accumulate(snapoffWettingLayerArea_.begin(), snapoffWettingLayerArea_.end(), 0.0), totalThroatCrossSectionalArea);
        regularBoundaryWettingThroatAreaSnapoff_[1] = std::min(std::accumulate(deltaPcSnapoffWettingLayerArea_.begin(), deltaPcSnapoffWettingLayerArea_.end(), 0.0), totalThroatCrossSectionalArea);
        regularBoundaryNonWettingThroatAreaSnapoff_[0] = totalThroatCrossSectionalArea - regularBoundaryWettingThroatAreaSnapoff_[0];
        regularBoundaryNonWettingThroatAreaSnapoff_[1] = totalThroatCrossSectionalArea - regularBoundaryWettingThroatAreaSnapoff_[1];

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

    // two parameters to specify regularization
    // regularization interval shows the range/width
    // regularization position shows the percetage of interval above critical pressure
    // e.g. for invasion event with entry pressure pce
    // the regularization area is [(pce+(position-1)*interval), (pce+position*interval)]
    /*!
     * \brief Returns the regularization Interval for invasion event
     */
    Scalar regularizationPcEntry(const int Idx) const
    { return regularizationPcEntry_[Idx]; }

    /*!
     * \brief Returns the regularization Interval for snap-off event
     */
    Scalar regularizationPcSnapoff(const int Idx) const
    { return regularizationPcSnapoff_[Idx]; }

    /*!
     * \brief Returns the curvature radius within the throat corresponds to entry pressure.
     */
    Scalar curvatureRadiusEntry(const int Idx) const
    { return curvatureRadiusEntry_[Idx];}

    /*!
     * \brief Returns the curvature radius within the throat corresponds to snap off pressure.
     */
    Scalar curvatureRadiusSnapoff(const int Idx) const
    { return curvatureRadiusSnapoff_[Idx];}

    /*!
     * \brief Returns delta pc used to calculate numerical derivative of K
     */
    Scalar deltaPc() const
    { return deltaPc_; }

    /*!
     * \brief Returns wetting layer area for each corner
     * at regularization pc right interval boundary for invasion
     */
    Scalar entryWettingLayerArea(const int cornerIdx) const
    { return entryWettingLayerArea_[cornerIdx]; }

    Scalar deltaPcEntryWettingLayerArea(const int cornerIdx) const
    { return deltaPcEntryWettingLayerArea_[cornerIdx]; }

    /*!
     * \brief Returns wetting layer area for each corner
     * at regularization pc right interval boundary for snap-off
     */
    Scalar snapoffWettingLayerArea(const int cornerIdx) const
    { return snapoffWettingLayerArea_[cornerIdx]; }

    Scalar deltaPcSnapoffWettingLayerArea(const int cornerIdx) const
    { return deltaPcSnapoffWettingLayerArea_[cornerIdx]; }

    /*!
     * \brief Returns total wetting layer area for invasion
     */
    Scalar regularBoundaryWettingThroatAreaEntry(const int Idx) const
    { return regularBoundaryWettingThroatAreaEntry_[Idx]; }

    /*!
     * \brief Returns total wetting layer area for snap-off
     */
    Scalar regularBoundaryWettingThroatAreaSnapoff(const int Idx) const
    { return regularBoundaryWettingThroatAreaSnapoff_[Idx]; }

    /*!
     * \brief Returns total non-wetting throat area for invasion
     */
    Scalar regularBoundaryNonWettingThroatAreaEntry(const int Idx) const
    { return regularBoundaryNonWettingThroatAreaEntry_[Idx]; }

    /*!
     * \brief Returns total non-wetting throat area for snap-off
     */
    Scalar regularBoundaryNonWettingThroatAreaSnapoff(const int Idx) const
    { return regularBoundaryNonWettingThroatAreaSnapoff_[Idx]; }

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

    // regularization parameters
    Scalar deltaPc_; // the delta pc used to calculate derivatives
    Scalar regularPcInterval_, regularizationPosition_; // regularization width and position
    Scalar regSaturationPercentage_;
    std::array<Scalar, 3> regularizationPcEntry_; // left and right Pc at interval for invasion
    std::array<Scalar, 3> regularizationPcSnapoff_; // left and right Pc at interval for snap-off
    std::array<Scalar, 2> curvatureRadiusEntry_; // curvature for entry
    std::array<Scalar, 2> curvatureRadiusSnapoff_; // curvature for snap-off
    NumCornerVector entryWettingLayerArea_, deltaPcEntryWettingLayerArea_;
    NumCornerVector snapoffWettingLayerArea_, deltaPcSnapoffWettingLayerArea_;
    std::array<Scalar, 2> regularBoundaryWettingThroatAreaEntry_;
    std::array<Scalar, 2> regularBoundaryWettingThroatAreaSnapoff_;
    std::array<Scalar, 2> regularBoundaryNonWettingThroatAreaEntry_;
    std::array<Scalar, 2> regularBoundaryNonWettingThroatAreaSnapoff_;
    bool regularizeWithSaturation_;
};


template<class AdvectionType>
class TwoPFluxVariablesCache<AdvectionType, 4, false>
{
    using Scalar = typename AdvectionType::Scalar;
    static constexpr auto numPhases = 2;
    using NumCornerVector = Dune::ReservedVector<Scalar, 4>;

public:

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

        if (invaded) // two-phase flow
        {
            const Scalar theta = spatialParams.contactAngle(element, elemVolVars);
            for (int i = 0; i< cornerHalfAngles.size(); ++i)
                wettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(curvatureRadius(), theta, cornerHalfAngles[i]);

            // make sure the wetting phase area does not exceed the total cross-section area
            throatCrossSectionalArea_[wPhaseIdx()] = std::min(
                std::accumulate(wettingLayerArea_.begin(), wettingLayerArea_.end(), 0.0),
                totalThroatCrossSectionalArea
            );
            throatCrossSectionalArea_[nPhaseIdx()] = totalThroatCrossSectionalArea - throatCrossSectionalArea_[wPhaseIdx()];
        }
        else // single-phase flow
        {
            for (int i = 0; i< cornerHalfAngles.size(); ++i)
                wettingLayerArea_[i] = 0.0;

            throatCrossSectionalArea_[wPhaseIdx()] = totalThroatCrossSectionalArea;
            throatCrossSectionalArea_[nPhaseIdx()] = 0.0;
        }

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
