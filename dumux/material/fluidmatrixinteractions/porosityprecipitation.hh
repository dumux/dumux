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
 * \ingroup Fluidmatrixinteractions
 * \brief Class for the evaluation of the porosity subject to precipitation.
 */
#ifndef DUMUX_POROSITY_PRECIPITATION_HH
#define DUMUX_POROSITY_PRECIPITATION_HH

#include <dumux/discretization/evalsolution.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Calculates the porosity depending on the volume fractions of precipitated minerals.
 *
 * \tparam Scalar The type used for scalar values
 * \param numComp The number of components in the fluid phases
 * \param numSolidPhases The number of precipitating solid phases
 */
template<class Scalar, int numComp, int numSolidPhases>
class PorosityPrecipitation
{
public:
    /*!
     * \brief Calculates the porosity in a sub-control volume
     * \param element Element
     * \param elemSol The element solution
     * \param scv Sub control volume
     * \param refPoro The solid matrix porosity without precipitates
     * \param minPoro A minimum porosity value
     */
    template<class Element, class SubControlVolume, class ElemSol>
    Scalar evaluatePorosity(const Element& element,
                            const SubControlVolume& scv,
                            const ElemSol& elemSol,
                            Scalar refPoro,
                            Scalar minPoro = 0.0) const
    {
        auto priVars = evalSolution(element, element.geometry(), elemSol, scv.center());

        Scalar sumPrecipitates = 0.0;
        for (unsigned int solidPhaseIdx = 0; solidPhaseIdx < numSolidPhases; ++solidPhaseIdx)
            sumPrecipitates += priVars[numComp + solidPhaseIdx];

        using std::max;
        return max(minPoro, refPoro - sumPrecipitates);
    }
};

} // namespace Dumux

#endif
