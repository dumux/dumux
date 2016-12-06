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
 * \brief  Defines the indices for the one-phase fully implicit model.
 */
#ifndef DUMUX_1P_INDICES_HH
#define DUMUX_1P_INDICES_HH

namespace Dumux
{
// \{
/*!
 * \ingroup NavierStokesModel
 * \ingroup ImplicitIndices
 * \brief The common indices for the isothermal stokes model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset = 0>
struct NavierStokesCommonIndices
{
    // number of dimensions
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    static const int dim = Grid::dimension;

    // Primary variable indices
    static const int momentumXIdx = PVOffset + 0; //!< Index of the x-component of the momentum equation
    static const int momentumYIdx = PVOffset + 1; //!< Index of the y-component of the momentum equation
    static const int momentumZIdx = PVOffset + 2; //!< Index of the z-component of the momentum equation
    static const int lastMomentumIdx = momentumXIdx+dim-1; //!< Index of the last component of the momentum equation

    static const int dimXIdx = 0; //!< Index of the x-component of a vector of size dim
    static const int dimYIdx = 1; //!< Index of the y-component of a vector of size dim
    static const int dimZIdx = 2; //!< Index of the z-component of a vector of size dim

    static const int massBalanceIdx = PVOffset + 0; //!< Index of the mass balance equation
    static const int conti0EqIdx = massBalanceIdx; //!< Index of first (for C-guys: 0th) mass conservation equation

    static const int pressureIdx = massBalanceIdx; //!< Index of the pressure in a solution vector
    static const int velocityXIdx = momentumXIdx; //!< Index of the x-component of the velocity
    static const int velocityYIdx = momentumYIdx; //!< Index of the y-component of the velocity
    static const int velocityZIdx = momentumZIdx; //!< Index of the z-component of the velocity

    static const int phaseIdx = 0; //!< Index of the fluid phase (required to use the same fluid system in coupled models)
};

// \}
} // end namespace

#endif
