// $Id: richardsproperties.hh 3840 2010-07-15 10:14:15Z bernd $
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief Index names for the Richards model.
 */
#ifndef DUMUX_RICHARDS_INDICES_HH
#define DUMUX_RICHARDS_INDICES_HH

namespace Dumux
{
/*!
 * \addtogroup RichardsModel
 */
// \{

/*!
 * \brief Index names for the Richards model.
 */
struct RichardsIndices
{
    //////////
    // primary variable indices
    //////////

    //! Primary variable index for the wetting phase pressure
    static const int pwIdx = 0;

    //////////
    // equation indices
    //////////
    //! Equation index for the mass conservation of the wetting phase
    static const int contiEqIdx = 0;

    //////////
    // phase indices
    //////////
    //! Phase index for the wetting phase
    static const int wPhaseIdx = 0;
    //! Phase index for the wetting phase
    static const int nPhaseIdx = 1;
};
// \}

} // end namepace

#endif
