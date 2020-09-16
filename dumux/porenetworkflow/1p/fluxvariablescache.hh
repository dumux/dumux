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
 * \brief Flux variables cache for 1p PNM
 */
#ifndef DUMUX_PNM_1P_FLUXVARIABLESCACHE_HH
#define DUMUX_PNM_1P_FLUXVARIABLESCACHE_HH

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
template<class AdvectionType>
class PNMOnePFluxVariablesCache
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
        throatRadius_ = problem.spatialParams().throatRadius(element, elemVolVars);

        cache_.fill(problem, element, fvGeometry, scvf, elemVolVars, *this, /*phaseIdx*/0);
        transmissibility_ = AdvectionType::calculateTransmissibility(problem, element, fvGeometry, scvf, elemVolVars, *this, /*phaseIdx*/0);
    }

    Throat::Shape throatCrossSectionShape() const
    { return throatCrossSectionShape_; }

    Scalar throatShapeFactor() const
    { return throatShapeFactor_; }

    Scalar transmissibility(const int phaseIdx = 0) const
    { return transmissibility_; }

    Scalar throatCrossSectionalArea(const int phaseIdx = 0) const
    { return throatCrossSectionalArea_; }

    Scalar throatLength() const
    { return throatLength_; }

    Scalar throatRadius() const
    { return throatRadius_; }

    const auto& singlePhaseFlowVariables() const
    { return cache_; }

private:
    Throat::Shape throatCrossSectionShape_;
    Scalar throatShapeFactor_;
    Scalar transmissibility_;
    Scalar throatCrossSectionalArea_;
    Scalar throatLength_;
    Scalar throatRadius_;

    typename AdvectionType::Transmissibility::SinglePhaseCache cache_;
};

} // end namespace

#endif
