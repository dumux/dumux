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
 * \ingroup SolidStates
 * \brief Update the solid volume fractions (inert and reacitve) and set them in the solidstate
 */
#ifndef DUMUX_UPDATE_SOLID_VOLUME_FRACTION_HH
#define DUMUX_UPDATE_SOLID_VOLUME_FRACTION_HH

#include <dune/common/concept.hh>
#include <dumux/discretization/solutionstate.hh>

namespace Dumux {

/*!
 * \ingroup SolidStates
 * \brief update the solid volume fractions (inert and reacitve) and set them in the solidstate
 * \note updates the inert components (TODO: these are assumed to come last right now in the solid system!)
 * \note gets the non-inert components from the primary variables
 */
template<class ElemState, class ExtVariables, class Problem, class Element, class Scv, class SolidState>
void updateSolidVolumeFractions(const ElemState& elemState,
                                const ExtVariables& extVariables,
                                const Problem& problem,
                                const Element& element,
                                const Scv& scv,
                                SolidState& solidState,
                                const int solidVolFracOffset)
{
    // compatibility layer with new elemsol-state + external variables style
    static constexpr bool isElemState = Dune::models<Experimental::Concept::ElementSolutionState, ElemState>();

    const auto& priVars = [&elemState, &scv] ()
    {
        if constexpr (isElemState)
            return elemState.elementSolution()[scv.localDofIndex()];
        else
            return elemState[scv.localDofIndex()];
    } ();

    for (int sCompIdx = solidState.numComponents-solidState.numInertComponents; sCompIdx < solidState.numComponents; ++sCompIdx)
    {
        const auto& sp = problem.spatialParams();
        using SolidSystem = typename SolidState::SolidSystem;
        if constexpr (isElemState)
        {
            const auto inertVolumeFraction = sp.template inertVolumeFraction<SolidSystem>(element, scv, elemState, extVariables, sCompIdx);
            solidState.setVolumeFraction(sCompIdx, inertVolumeFraction);
        }
        else
        {
            const auto inertVolumeFraction = sp.template inertVolumeFraction<SolidSystem>(element, scv, elemState, sCompIdx);
            solidState.setVolumeFraction(sCompIdx, inertVolumeFraction);
        }
    }

    if (!(solidState.isInert()))
    {
        for (int sCompIdx = 0; sCompIdx < solidState.numComponents-solidState.numInertComponents; ++sCompIdx)
            solidState.setVolumeFraction(sCompIdx, priVars[solidVolFracOffset + sCompIdx]);
    }
}

/*!
 * \ingroup SolidStates
 * \brief update the solid volume fractions (inert and reacitve) and set them in the solidstate
 * \note updates the inert components (TODO: these are assumed to come last right now in the solid system!)
 * \note gets the non-inert components from the primary variables
 */
template<class ElemState, class Problem, class Element, class Scv, class SolidState>
void updateSolidVolumeFractions(const ElemState& elemState,
                                const Problem& problem,
                                const Element& element,
                                const Scv& scv,
                                SolidState& solidState,
                                const int solidVolFracOffset)
{
    struct EmptyExtVars {} empty;
    updateSolidVolumeFractions(elemState, empty, problem, element, scv, solidState, solidVolFracOffset);
}

} // end namespace Dumux

#endif
