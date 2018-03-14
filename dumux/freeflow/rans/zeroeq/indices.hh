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
 * \ingroup ZeroEqModel
 * \copydoc Dumux::ZeroEqIndices
 */
#ifndef DUMUX_ZEROEQ_INDICES_HH
#define DUMUX_ZEROEQ_INDICES_HH

#include <dumux/freeflow/navierstokes/indices.hh>

namespace Dumux
{
// \{
/*!
 * \ingroup ZeroEqModel
 * \brief The common indices for the isothermal ZeroEq model.
 *
 * \tparam dimension The dimension of the problem
 * \tparam numEquations the number of model equations
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <int dimension, int numEquations, int PVOffset = 0>
struct ZeroEqIndices
    : NavierStokesIndices<dimension, numEquations, PVOffset>
{
    static constexpr int noEddyViscosityModel = 0;
    static constexpr int prandtl = 1;
    static constexpr int modifiedVanDriest = 2;
};

// \}
} // end namespace

#endif
