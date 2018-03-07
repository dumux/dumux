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
 * \ingroup NavierStokesNIModel
 * \copydoc Dumux::NavierStokesNonIsothermalIndices
 */
#ifndef DUMUX_NAVIERSTOKES_NI_INDICES_HH
#define DUMUX_NAVIERSTOKES_NI_INDICES_HH

#include <dumux/common/properties.hh>

namespace Dumux
{
// \{
/*!
 * \ingroup NavierStokesNIModel
 * \brief Indices for the non-isothermal Navier-Stokes model.
 *
 * \tparam dimension The dimension of the problem
 * \tparam numEquations The number of model equations
 * \tparam IsothermalIndices The isothermal indices class
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <int dimension, int numEquations, class IsothermalIndices, int PVOffset = 0>
class NavierStokesNonIsothermalIndices : public IsothermalIndices
{
public:
    static constexpr auto dim = dimension;
    static constexpr auto numEq = numEquations;

    static constexpr auto energyBalanceIdx = PVOffset + numEq - dim - 1;
    static constexpr int temperatureIdx = energyBalanceIdx;
};
} // end namespace

#endif // DUMUX_NAVIERSTOKES_NI_INDICES_HH
