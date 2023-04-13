// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMOnePModel
 * \copydoc Dumux::PoreNetwork::OnePFluxVariablesCache
 */
#ifndef DUMUX_PNM_1P_FLUXVARIABLESCACHE_HH
#define DUMUX_PNM_1P_FLUXVARIABLESCACHE_HH

#include <dumux/porenetwork/common/throatproperties.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMOnePModel
 * \brief Flux variables cache for the single-phase-flow PNM
 *        Store data required for flux calculation.
 */
template<class AdvectionType>
class OnePFluxVariablesCache
{
    using Scalar = typename AdvectionType::Scalar;
public:
    //! whether the cache needs an update when the solution changes
    static bool constexpr isSolDependent = true;

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
