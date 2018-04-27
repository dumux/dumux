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
 * \ingroup SolidStates
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase fluid system assuming immiscibility and
 *        thermodynamic equilibrium.
 */
#ifndef DUMUX_BASE_SOLID_STATE_HH
#define DUMUX_BASE_SOLID_STATE_HH

#include <dumux/common/valgrind.hh>

#include <limits>

namespace Dumux
{

template<class ElemSol, class Problem, class Element, class Scv, class SolidState>
void updateSolidVolumeFractions(const ElemSol &elemSol,
                                const Problem &problem,
                                const Element &element,
                                const Scv &scv,
                                SolidState& solidState,
                                const int numFluidComponents)
{
    for(int sCompIdx = solidState.numComponents- solidState.numInertComponents; sCompIdx < solidState.numComponents; ++sCompIdx)
    {
        solidState.setVolumeFraction(sCompIdx, problem.spatialParams().inertVolumeFraction(element, scv, solidState, sCompIdx));
    }
    if (!(solidState.isInert()))
    {
        auto&& priVars = elemSol[scv.localDofIndex()];
        for(int sCompIdx = 0; sCompIdx < solidState.numComponents- solidState.numInertComponents; ++sCompIdx)
        {
           solidState.setVolumeFraction(sCompIdx, priVars[numFluidComponents + sCompIdx]);
        }
    }
}

} // end namespace Dumux

#endif
