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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesIndices
 */
#ifndef DUMUX_NAVIERSTOKES_INDICES_HH
#define DUMUX_NAVIERSTOKES_INDICES_HH

#include <dumux/common/properties.hh>

namespace Dumux {

// \{
/*!
 * \ingroup NavierStokesModel
 * \brief The common indices for the isothermal Navier-Stokes model.
 *
 * \tparam dimension The dimension of the problem
 * \tparam numEquations the number of model equations
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <int dimension, int numEquations, int PVOffset = 0>
struct NavierStokesIndices
{

    static constexpr int dimXIdx = 0; //!< Index of the x-component of a vector of size dim
    static constexpr int dimYIdx = 1; //!< Index of the y-component of a vector of size dim
    static constexpr int dimZIdx = 2; //!< Index of the z-component of a vector of size dim

    static constexpr int massBalanceIdx = PVOffset + 0; //!< Index of the mass balance equation
    static constexpr int conti0EqIdx = massBalanceIdx; //!< Index of the mass balance equation
    static constexpr int pressureIdx = massBalanceIdx; //!< Index of the pressure in a solution vector

    static constexpr auto dim = dimension;
    static constexpr auto numEq = numEquations;
    static constexpr auto momentumBalanceOffset = numEquations - dim;

    static constexpr int momentumBalanceIdx = PVOffset + momentumBalanceOffset; //!< Index of the momentum balance equation
    static constexpr int momentumXBalanceIdx = momentumBalanceIdx; //!< Index of the momentum balance equation
    static constexpr int momentumYBalanceIdx = momentumBalanceIdx + 1; //!< Index of the momentum balance equation
    static constexpr int momentumZBalanceIdx = momentumBalanceIdx + 2; //!< Index of the momentum balance equation

    static constexpr int velocityXIdx = momentumBalanceIdx; //!< Index of the velocity in a solution vector
    static constexpr int velocityYIdx = momentumBalanceIdx + 1; //!< Index of the velocity in a solution vector
    static constexpr int velocityZIdx = momentumBalanceIdx + 2; //!< Index of the velocity in a solution vector

    /*!
     * \brief Index of the velocity in a solution vector given a certain direction.
     *
     * \param dirIdx The index of the direction.
     */
    static constexpr int velocity(int dirIdx)
    {
        return dirIdx + momentumBalanceIdx;
    }
};

// \}
} // end namespace

#endif
