// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
 /*!
  * \file
  * \ingroup PNMTwoPModel
  * \brief Flux variables cache for the two-phase-flow PNM with a regularized throat state
  */
#ifndef DUMUX_PNM_2P_REGULARIZED_FLUXVARIABLESCACHE_HH
#define DUMUX_PNM_2P_REGULARIZED_FLUXVARIABLESCACHE_HH

#include <array>
#include <algorithm>
#include <numeric>
#include <dune/common/reservedvector.hh>
#include <dune/common/std/type_traits.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>
#include <dumux/porenetwork/common/throatproperties.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMTwoPModel
 * \brief Flux variables cache for the two-phase-flow PNM with a regularized throat state
 *
 * In addition to the quantities of TwoPFluxVariablesCache, this cache carries the throat state
 * theta, the fraction of the time step at which the invasion or the snap-off of a throat occurs.
 * Theta is a regularized Heaviside function of the capillary pressure, around the entry pressure
 * if the throat is not invaded and around the snap-off pressure if it is. The width of the
 * regularization is set by "Problem.RegularizationDelta".
 *
 * The throat transmissibilities are weighted with theta rather than being switched discretely, so
 * the wetting layer areas at the entry pressure are cached as well. The discrete invasion state is
 * only updated once the time step has converged, see StateSwitchMethod::EndOfTimeStep.
 */
template<class AdvectionType, int maxNumCorners = 4>
class TwoPRegularizedFluxVariablesCache
{
    using ThisType = TwoPRegularizedFluxVariablesCache<AdvectionType, maxNumCorners>;
    using Scalar = typename AdvectionType::Scalar;
    static constexpr auto numPhases = 2;
    using NumCornerVector = Dune::ReservedVector<Scalar, maxNumCorners>;

public:
    //! whether the cache needs an update when the solution changes
    static bool constexpr isSolDependent = true;

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
        const Scalar totalThroatCrossSectionalArea = spatialParams.throatCrossSectionalArea(element, elemVolVars);
        const Scalar alpha = spatialParams.contactAngle(element, elemVolVars);

        wettingLayerArea_.clear(); wettingLayerArea_.resize(cornerHalfAngles.size());
        entryWettingLayerArea_.clear(); entryWettingLayerArea_.resize(cornerHalfAngles.size());

        for (int i = 0; i < cornerHalfAngles.size(); ++i)
        {
            wettingLayerArea_[i] = std::min(
                Throat::wettingLayerCrossSectionalArea(curvatureRadius(), alpha, cornerHalfAngles[i]),
                totalThroatCrossSectionalArea
            );
            entryWettingLayerArea_[i] = Throat::wettingLayerCrossSectionalArea(curvatureRadiusInvasion(), alpha, cornerHalfAngles[i]);
        }

        const auto wettingAreaEntry = std::min(
            std::accumulate(entryWettingLayerArea_.begin(), entryWettingLayerArea_.end(), 0.0),
            totalThroatCrossSectionalArea
        );
        nonWettingAreaEntry_ = totalThroatCrossSectionalArea - wettingAreaEntry;

        // make sure the wetting phase area does not exceed the total cross-section area
        throatCrossSectionalArea_[wPhaseIdx()] = std::min(
            std::accumulate(wettingLayerArea_.begin(), wettingLayerArea_.end(), 0.0),
            totalThroatCrossSectionalArea
        );
        throatCrossSectionalArea_[nPhaseIdx()] = totalThroatCrossSectionalArea - throatCrossSectionalArea_[wPhaseIdx()];

        // the throat state is evaluated once the quantities above are available
        static const Scalar regularizationDelta = getParamFromGroup<Scalar>(problem.paramGroup(), "Problem.RegularizationDelta", 1e-2);
        theta_ = computeTheta_(regularizationDelta);

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
    { return surfaceTension_ / std::max(pc_, 1e-16); }

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
     * \brief Returns the throat state, i.e. the fraction of the time step at which
     *        the invasion or the snap-off of the throat occurs, within [0,1].
     */
    Scalar theta() const
    { return theta_; }

    //! Returns the curvature radius within the throat at the entry pressure
    Scalar curvatureRadiusInvasion() const
    { return surfaceTension_/pcEntry_; }

    //! Returns the wetting layer area of a corner at the entry pressure
    Scalar entryWettingLayerArea(const int cornerIdx) const
    { return entryWettingLayerArea_[cornerIdx]; }

    //! Returns the nonwetting phase cross-sectional area at the entry pressure
    Scalar nonWettingAreaEntry() const
    { return nonWettingAreaEntry_; }

private:

    /*!
     * \brief Returns the regularized throat state theta, see the class documentation.
     */
    Scalar computeTheta_(const Scalar delta) const
    {
        using std::abs; using std::sin; using std::clamp;

        if (!invaded_)
        {
            // regularize once pc is above the entry pressure
            const Scalar dp = clamp(pc_/pcEntry_ - Scalar(1.0), Scalar(0.0), delta);
            return Scalar(0.5*(1.0 + sin(M_PI*(dp - 0.5*delta)/delta)));
        }
        else
        {
            // regularize once pc is below the snap-off pressure
            // a vanishing snap-off pressure means that snap-off cannot occur
            if (pcSnapoff_ == Scalar(0.0))
                return Scalar(1.0);

            const Scalar dp = clamp(pc_/abs(pcSnapoff_) - sign(pcSnapoff_), -delta, Scalar(0.0));
            return Scalar(0.5*(1.0 + sin(M_PI*(dp + 0.5*delta)/delta)));
        }
    }

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
    NumCornerVector entryWettingLayerArea_;
    Scalar nonWettingAreaEntry_;
    Scalar theta_;
    std::size_t nPhaseIdx_;
    Scalar poreToPoreDistance_;

    typename AdvectionType::Transmissibility::SinglePhaseCache singlePhaseCache_;
    typename AdvectionType::Transmissibility::NonWettingPhaseCache nonWettingPhaseCache_;
    typename AdvectionType::Transmissibility::WettingLayerCache wettingLayerCache_;
};

} // end Dumux::PoreNetwork

#endif
