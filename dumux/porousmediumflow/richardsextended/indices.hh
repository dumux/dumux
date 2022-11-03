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
 * \ingroup ExtendedRichardsModel
 * \brief Index names for the extended Richards model.
 */

#ifndef DUMUX_RICHARDSEXTENDED_INDICES_HH
#define DUMUX_RICHARDSEXTENDED_INDICES_HH

#include <dumux/porousmediumflow/richards/indices.hh>

namespace Dumux {

/*!
 * \ingroup ExtendedRichardsModel
 * \brief Index names for the extended Richards model.
 */

struct ExtendedRichardsIndices : public RichardsIndices
{
    static constexpr int switchIdx = 0;

    // present phases (-> 'pseudo' primary variable)
    static constexpr int liquidPhaseOnly = 1; //!< Only the liquid phase is present
    static constexpr int gasPhaseOnly = 2; //!< Only the gas phase is present
    static constexpr int bothPhases = 3; //!< Both phases are present
};

} // end namespace Dumux

#endif
