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
 * \ingroup RichardsNCModel
 * \brief Defines the primary variable and equation indices used by
 *        the richardsnc model.
 */

#ifndef DUMUX_RICHARDSNC_INDICES_HH
#define DUMUX_RICHARDSNC_INDICES_HH

namespace Dumux
{
// \{

/*!
 * \ingroup RichardsNCModel
 * \brief The indices for the isothermal Richards, n-component model.
 */
template <int PVOffset = 0>
struct RichardsNCIndices
{

    //! Set the index of the phases for accessing the volvars
    static const int wPhaseIdx = 0;
    static const int nPhaseIdx = 1;

    //! Component indices
    static const int compMainIdx = PVOffset + 0; //!< main component index

    //! primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< pressure

    //! \note These indices make sense if the first balance is replaced by the
    //!       total mass balance.

    //! Equation indices
    static const int conti0EqIdx = PVOffset + 0; //!< continuity equation index
};

// \}
}

#endif
