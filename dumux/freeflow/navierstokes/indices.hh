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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesIndices
 */
#ifndef DUMUX_NAVIERSTOKES_INDICES_HH
#define DUMUX_NAVIERSTOKES_INDICES_HH

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The common indices for the isothermal Navier-Stokes model.
 *
 * \tparam dimension The dimension of the problem
 */
template <int dimension>
struct NavierStokesIndices
{
    static constexpr int dimXIdx = 0; //!< Index of the x-component of a vector of size dim
    static constexpr int dimYIdx = 1; //!< Index of the y-component of a vector of size dim
    static constexpr int dimZIdx = 2; //!< Index of the z-component of a vector of size dim

    static constexpr int conti0EqIdx = dimension; //!< Index of the total mass balance equation
    static constexpr int pressureIdx = conti0EqIdx; //!< Index of the pressure in a solution vector

    static constexpr auto dim = dimension;

    static constexpr int momentumXBalanceIdx = 0; //!< Index of the momentum balance equation
    static constexpr int momentumYBalanceIdx = 1; //!< Index of the momentum balance equation
    static constexpr int momentumZBalanceIdx = 2; //!< Index of the momentum balance equation

    static constexpr int velocityXIdx = 0; //!< Index of the velocity in a solution vector
    static constexpr int velocityYIdx = 1; //!< Index of the velocity in a solution vector
    static constexpr int velocityZIdx = 2; //!< Index of the velocity in a solution vector

    /*!
     * \brief Index of the velocity in a solution vector given a certain direction.
     *
     * \param dirIdx The index of the direction.
     */
    static constexpr int velocity(int dirIdx)
    {
        return dirIdx;
    }
};

} // end namespace Dumux

#endif
