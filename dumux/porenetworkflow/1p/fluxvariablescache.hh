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
 * \brief Flux variables cache for the single-phase-flow PNM
 */
#ifndef DUMUX_PNM_1P_FLUXVARIABLESCACHE_HH
#define DUMUX_PNM_1P_FLUXVARIABLESCACHE_HH

#include <dumux/porenetworkflow/common/throatproperties.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PoreNetworkOnePModel
 * \brief Flux variables cache for the single-phase-flow PNM
 *        Store data required for flux calculation.
 */
template<class AdvectionType>
class OnePFluxVariablesCache
{
    using Scalar = typename AdvectionType::Scalar;
public:

    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables>
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const typename FVElementGeometry::SubControlVolumeFace& scvf)
    {
        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
        throatCrossSectionShape_ = fvGeometry.gridGeometry().throatCrossSectionShape(eIdx);
        throatShapeFactor_ = fvGeometry.gridGeometry().throatShapeFactor(eIdx);
        throatCrossSectionalArea_ = problem.spatialParams().throatCrossSectionalArea(element, elemVolVars);
        throatLength_ = problem.spatialParams().throatLength(element, elemVolVars);
        throatInscribedRadius_ = problem.spatialParams().throatInscribedRadius(element, elemVolVars);
        poreToPoreDistance_ = element.geometry().volume();

        cache_.fill(problem, element, fvGeometry, scvf, elemVolVars, *this, /*phaseIdx*/0);
        transmissibility_ = AdvectionType::calculateTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, *this, /*phaseIdx*/0);
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
    Scalar transmissibility(const int phaseIdx = 0) const
    { return transmissibility_; }

    /*!
     * \brief Returns the throats's cross-sectional area.
     */
    Scalar throatCrossSectionalArea(const int phaseIdx = 0) const
    { return throatCrossSectionalArea_; }

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
     * \brief Returns the throats's pore-to-pore-center distance.
     */
    Scalar poreToPoreDistance() const
    { return poreToPoreDistance_; }

    /*!
     * \brief Returns the throats's cached flow variables for single-phase flow.
     */
    const auto& singlePhaseFlowVariables() const
    { return cache_; }

private:
    Throat::Shape throatCrossSectionShape_;
    Scalar throatShapeFactor_;
    Scalar transmissibility_;
    Scalar throatCrossSectionalArea_;
    Scalar throatLength_;
    Scalar throatInscribedRadius_;
    Scalar poreToPoreDistance_;

    typename AdvectionType::Transmissibility::SinglePhaseCache cache_;
};

} // end namespace Dumux::PoreNetwork

#endif
