// $Id: 1p2cproperties.hh 3838 2010-07-15 08:31:53Z bernd $
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Defines the primary variable and equation indices used by
 *        the 1p2c model
 */

#ifndef DUMUX_1P2C_INDICES_HH
#define DUMUX_1P2C_INDICES_HH

namespace Dumux
{
/*!
 * \addtogroup OnePTwoCBoxModel
 */
// \{

/*!
 * \brief The indices for the isothermal single-phase, two-component model.
 */
template <int PVOffset = 0>
struct OnePTwoCIndices
{
    // Equation indices
    static const int contiEqIdx = PVOffset + 0; //!< continuity equation index
    static const int transEqIdx = PVOffset + 1; //!< transport equation index

    // primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< pressure
    static const int x1Idx = PVOffset + 1; //!< mole fraction of the second component
                                //in my case the therapeutic agent
};

// \}
}

#endif

