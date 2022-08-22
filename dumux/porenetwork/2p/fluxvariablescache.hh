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
        regPercent_ = getParam<Scalar>("Regularization.RegPercentage", 0.01);
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
        pcEntry_ = problem.spatialParams().pcEntry(element, elemVolVars);
        pcSnapoff_ = problem.spatialParams().pcSnapoff(element, elemVolVars);

        // take the average surface tension of both adjacent pores TODO: is this correct?
        surfaceTension_ = 0.5*(elemVolVars[0].surfaceTension() + elemVolVars[1].surfaceTension());

        pc_ = std::max(elemVolVars[0].capillaryPressure(), elemVolVars[1].capillaryPressure());

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& elemSol = elementSolution(element, elemVolVars, fvGeometry);
        const auto& spatialParams = problem.spatialParams();
        auto fluidMatrixInteraction = spatialParams.fluidMatrixInteraction(element, insideScv, elemSol);
        // open a file to write the analytical pcEntry and Sw for comparison purposes   
        std::ofstream file("analytical_SwEntry");
        file << std::setprecision(15) << "pcEntry: " << pcEntry_ << " Sw: " << fluidMatrixInteraction.sw(pcEntry_) << std::endl;

        regInvasionInterval_[0] = pcEntry_;
        regSnapoffInterval_[1] = pcSnapoff_;

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

        throatInscribedRadius_ = problem.spatialParams().throatInscribedRadius(element, elemVolVars);
        throatLength_ = problem.spatialParams().throatLength(element, elemVolVars);
        invaded_ = invaded;
        poreToPoreDistance_ = element.geometry().volume();

        // get the non-wetting phase index
        using FluidSystem = typename ElementVolumeVariables::VolumeVariables::FluidSystem;
        nPhaseIdx_ = 1 - spatialParams.template wettingPhase<FluidSystem>(element, elemVolVars);


        const auto& cornerHalfAngles = spatialParams.cornerHalfAngles(element);
        wettingLayerArea_.clear(); wettingLayerArea_.resize(cornerHalfAngles.size());
        const Scalar totalThroatCrossSectionalArea = spatialParams.throatCrossSectionalArea(element, elemVolVars);
        entryWettingLayerArea_.clear(); entryWettingLayerArea_.resize(cornerHalfAngles.size());
        snapoffWettingLayerArea_.clear(); snapoffWettingLayerArea_.resize(cornerHalfAngles.size());
        epsilonEntryWettingLayerArea_.clear(); epsilonEntryWettingLayerArea_.resize(cornerHalfAngles.size());
        epsilonSnapoffWettingLayerArea_.clear(); epsilonSnapoffWettingLayerArea_.resize(cornerHalfAngles.size());

        const Scalar theta = spatialParams.contactAngle(element, elemVolVars);
        for (int i = 0; i< cornerHalfAngles.size(); ++i)
        {
            wettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(curvatureRadius(), theta, cornerHalfAngles[i]);
            entryWettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(curvatureRadiusInvasion(0), theta, cornerHalfAngles[i]);
            snapoffWettingLayerArea_[i] =  Throat::wettingLayerCrossSectionalArea(curvatureRadiusSnapoff(0), theta, cornerHalfAngles[i]);
            epsilonEntryWettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(curvatureRadiusInvasion(1), theta, cornerHalfAngles[i]);
            epsilonSnapoffWettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(curvatureRadiusSnapoff(1), theta, cornerHalfAngles[i]);
        }

        // make sure the wetting phase area does not exceed the total cross-section area
        throatCrossSectionalArea_[wPhaseIdx()] = std::min(
            std::accumulate(wettingLayerArea_.begin(), wettingLayerArea_.end(), 0.0),
            totalThroatCrossSectionalArea
        );
        throatCrossSectionalArea_[nPhaseIdx()] = totalThroatCrossSectionalArea - throatCrossSectionalArea_[wPhaseIdx()];

        regBoundaryWettingThroatAreaInvasion_[0] = std::min(std::accumulate(entryWettingLayerArea_.begin(), entryWettingLayerArea_.end(), 0.0), totalThroatCrossSectionalArea);
        regBoundaryWettingThroatAreaInvasion_[1] = std::min(std::accumulate(epsilonEntryWettingLayerArea_.begin(), epsilonEntryWettingLayerArea_.end(), 0.0), totalThroatCrossSectionalArea);
        regBoundaryNonwettingThroatAreaInvasion_[0] = totalThroatCrossSectionalArea - regBoundaryWettingThroatAreaInvasion_[0];
        regBoundaryNonwettingThroatAreaInvasion_[1] = totalThroatCrossSectionalArea - regBoundaryWettingThroatAreaInvasion_[1];
        regBoundaryWettingThroatAreaSnapoff_[0] = std::min(std::accumulate(snapoffWettingLayerArea_.begin(), snapoffWettingLayerArea_.end(), 0.0), totalThroatCrossSectionalArea);
        regBoundaryWettingThroatAreaSnapoff_[1] = std::min(std::accumulate(epsilonSnapoffWettingLayerArea_.begin(), epsilonSnapoffWettingLayerArea_.end(), 0.0), totalThroatCrossSectionalArea);
        regBoundaryNonwettingThroatAreaSnapoff_[0] = totalThroatCrossSectionalArea - regBoundaryWettingThroatAreaSnapoff_[0];
        regBoundaryNonwettingThroatAreaSnapoff_[1] = totalThroatCrossSectionalArea - regBoundaryWettingThroatAreaSnapoff_[1];


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

    /*!
     * \brief Returns the regularization Interval for invasion event
     */
    Scalar regInvasionInterval(const int Idx) const
    { return regInvasionInterval_[Idx]; }

    /*!
     * \brief Returns the regularization Interval for snap-off event
     */
    Scalar regSnapoffInterval(const int Idx) const
    { return regSnapoffInterval_[Idx]; }

    /*!
     * \brief Returns the curvature radius within the throat corresponds to entry pressure.
     */
    Scalar curvatureRadiusInvasion(const int Idx) const
    { return curvatureRadiusInvasion_[Idx];}

    /*!
     * \brief Returns the curvature radius within the throat corresponds to snap off pressure.
     */
    Scalar curvatureRadiusSnapoff(const int Idx) const
    { return curvatureRadiusSnapoff_[Idx];}

    /*!
     * \brief Returns delta pc used to calculate numerical derivative of K
     */
    Scalar epsilonPc() const
    { return epsilonPc_; }

    /*!
     * \brief Returns wetting layer area for each corner
     * at regularization pc right interval boundary for invasion
     */
    Scalar entryWettingLayerArea(const int cornerIdx) const
    { return entryWettingLayerArea_[cornerIdx]; }

    Scalar epsilonEntryWettingLayerArea(const int cornerIdx) const
    { return epsilonEntryWettingLayerArea_[cornerIdx]; }

    /*!
     * \brief Returns wetting layer area for each corner
     * at regularization pc right interval boundary for snap-off
     */
    Scalar snapoffWettingLayerArea(const int cornerIdx) const
    { return snapoffWettingLayerArea_[cornerIdx]; }

    Scalar epsilonSnapoffWettingLayerArea(const int cornerIdx) const
    { return epsilonSnapoffWettingLayerArea_[cornerIdx]; }

    /*!
     * \brief Returns total wetting layer area for invasion
     */
    Scalar regBoundaryWettingThroatAreaInvasion(const int Idx) const
    { return regBoundaryWettingThroatAreaInvasion_[Idx]; }

    /*!
     * \brief Returns total wetting layer area for snap-off
     */
    Scalar regBoundaryWettingThroatAreaSnapoff(const int Idx) const
    { return regBoundaryWettingThroatAreaSnapoff_[Idx]; }

    /*!
     * \brief Returns total non-wetting throat area for invasion
     */
    Scalar regBoundaryNonwettingThroatAreaInvasion(const int Idx) const
    { return regBoundaryNonwettingThroatAreaInvasion_[Idx]; }

    /*!
     * \brief Returns total non-wetting throat area for snap-off
     */
    Scalar regBoundaryNonWettingThroatAreaSnapoff(const int Idx) const
    { return regBoundaryNonwettingThroatAreaSnapoff_[Idx]; }

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

    static constexpr double epsilonPc_ = 1e-8;
    Scalar regPercent_;
    std::array<Scalar, 3> regInvasionInterval_;
    std::array<Scalar, 3> regSnapoffInterval_;
    std::array<Scalar, 2> curvatureRadiusInvasion_;
    std::array<Scalar, 2> curvatureRadiusSnapoff_;
    NumCornerVector entryWettingLayerArea_, epsilonEntryWettingLayerArea_;
    NumCornerVector snapoffWettingLayerArea_, epsilonSnapoffWettingLayerArea_;
    std::array<Scalar, 2> regBoundaryWettingThroatAreaInvasion_;
    std::array<Scalar, 2> regBoundaryWettingThroatAreaSnapoff_;
    std::array<Scalar, 2> regBoundaryNonwettingThroatAreaInvasion_;
    std::array<Scalar, 2> regBoundaryNonwettingThroatAreaSnapoff_;
};

} // end Dumux::PoreNetwork

#endif
