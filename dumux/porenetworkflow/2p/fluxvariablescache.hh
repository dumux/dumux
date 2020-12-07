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
 * \brief Base class for the flux variables
 */
#ifndef DUMUX_PNM_2P_FLUXVARIABLESCACHE_HH
#define DUMUX_PNM_2P_FLUXVARIABLESCACHE_HH

#include <array>
#include <algorithm>
#include <dune/common/reservedvector.hh>
#include <dumux/porenetworkflow/common/throatproperties.hh>

namespace Dumux
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//! The cache is dependent on the active physical processes (advection, diffusion, heat conduction)
//! For each type of process there is a base cache storing the data required to compute the respective fluxes
//! Specializations of the overall cache are provided for combinations of processes
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*!
 * \ingroup ImplicitModel
 * \brief The flux variables cache classes for porous media.
 *        Store data required for flux calculation. For each type of physical process (advection, diffusion, heat conduction)
 *        there is a base cache storing the data required to compute the respective fluxes. Specializations of the overall
 *        cache class are provided for different combinations of processes.
 */

//! We only store discretization-related quantities for the box method.
//! Thus, we need no physics-dependent specialization.
template<class AdvectionType, int maxNumCorners = 4>
class PNMTwoPFluxVariablesCache
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
                bool invaded)
    {
        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
        throatCrossSectionShape_ = fvGeometry.gridGeometry().throatCrossSectionShape(eIdx);
        throatShapeFactor_ = fvGeometry.gridGeometry().throatShapeFactor(eIdx);
        pc_ = std::max(elemVolVars[0].capillaryPressure(), elemVolVars[1].capillaryPressure());
        pcEntry_ = problem.spatialParams().pcEntry(element, elemVolVars);
        pcSnapoff_ = problem.spatialParams().pcSnapoff(element, elemVolVars);
        throatRadius_ = problem.spatialParams().throatRadius(element, elemVolVars);
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
            transmissibility_[phaseIdx] = AdvectionType::calculateTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, *this, phaseIdx);
    }

    Throat::Shape throatCrossSectionShape() const
    { return throatCrossSectionShape_; }

    Scalar throatShapeFactor() const
    { return throatShapeFactor_; }

    Scalar transmissibility(const int phaseIdx) const
    { return transmissibility_[phaseIdx]; }

    Scalar throatCrossSectionalArea(const int phaseIdx) const
    { return throatCrossSectionalArea_[phaseIdx]; }

    Scalar throatCrossSectionalArea() const
    { return throatCrossSectionalArea_[0] + throatCrossSectionalArea_[1]; }

    Scalar throatLength() const
    { return throatLength_; }

    Scalar throatRadius() const
    { return throatRadius_; }

    Scalar pcEntry() const
    { return pcEntry_; }

    Scalar pcSnapoff() const
    { return pcSnapoff_; }

    Scalar pc() const
    { return pc_; }

    Scalar surfaceTension() const
    { return surfaceTension_; }

    bool invaded() const
    { return invaded_; }

    Scalar curvatureRadius() const
    { return surfaceTension_ / pc_;}

    Scalar wettingLayerCrossSectionalArea(const int cornerIdx) const
    { return wettingLayerArea_[cornerIdx]; }

    std::size_t wPhaseIdx() const
    { return 1 - nPhaseIdx_; }

    std::size_t nPhaseIdx() const
    { return nPhaseIdx_; }

    const auto& singlePhaseFlowVariables() const
    { return singlePhaseCache_; }

    const auto& nonWettingPhaseFlowVariables() const
    { return nonWettingPhaseCache_; }

    const auto& wettingLayerFlowVariables() const
    { return wettingLayerCache_; }

    Scalar poreToPoreDistance() const
    { return poreToPoreDistance_; }

private:
    Throat::Shape throatCrossSectionShape_;
    Scalar throatShapeFactor_;
    std::array<Scalar, numPhases> transmissibility_;
    std::array<Scalar, numPhases> throatCrossSectionalArea_;
    Scalar throatLength_;
    Scalar throatRadius_;
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

} // end namespace

#endif
