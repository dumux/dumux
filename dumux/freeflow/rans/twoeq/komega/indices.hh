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
 * \ingroup KOmegaModel
 * \copydoc Dumux::KOmegaIndices
 */
#ifndef DUMUX_KOMEGA_INDICES_HH
#define DUMUX_KOMEGA_INDICES_HH

#include <dumux/freeflow/navierstokes/indices.hh>

namespace Dumux {

// \{
/*!
 * \ingroup KOmegaModel
 * \brief The common indices for the isothermal KOmega model.
 *
 * \tparam dimension The dimension of the problem
 * \tparam numComponents The number of considered transported components
 */
template<int dimension, int numComponents>
struct KOmegaIndices : public NavierStokesIndices<dimension>
{
public:
    static constexpr auto turbulentKineticEnergyEqIdx = dimension + numComponents;
    static constexpr auto turbulentKineticEnergyIdx = turbulentKineticEnergyEqIdx;
    static constexpr auto dissipationEqIdx = turbulentKineticEnergyEqIdx + 1;
    static constexpr auto dissipationIdx = dissipationEqIdx;
};

// \}
} // end namespace

#endif
