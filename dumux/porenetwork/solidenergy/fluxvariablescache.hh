// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup PNMSolidEnergyModel
 * \brief Flux variables cache for PNM solid-energy model
 */
#ifndef DUMUX_PNM_SOLID_ENERGY_FLUXVARIABLESCACHE_HH
#define DUMUX_PNM_SOLID_ENERGY_FLUXVARIABLESCACHE_HH

namespace Dumux::PoreNetwork {

template<class Scalar>
class SolidEnergyFluxVariablesCache
{
public:
    //! whether the cache needs an update when the solution changes
    static bool constexpr isSolDependent = false;

    template<class Problem, class Element, class FVElementGeometry,
             class ElementVolumeVariables, class SubControlVolumeFace>
    void update(const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace& scvf)
    {
        grainContactArea_ = problem.spatialParams().throatCrossSectionalArea(element, elemVolVars);
        throatLength_ = problem.spatialParams().throatLength(element, elemVolVars);
        throatInscribedRadius_ = problem.spatialParams().throatInscribedRadius(element, elemVolVars);
    }

    Scalar grainContactArea() const
    { return grainContactArea_; }

    Scalar throatLength() const
    { return throatLength_; }

    Scalar throatInscribedRadius() const
    { return throatInscribedRadius_; }

private:
    Scalar grainContactArea_;
    Scalar throatLength_;
    Scalar throatInscribedRadius_;
};

} // end namespace Dumux::PoreNetwork

#endif
