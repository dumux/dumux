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
 * \brief  Defines labels for pores and throats.
 */
#ifndef DUMUX_PNM_LABELS_HH
#define DUMUX_PNM_LABELS_HH

namespace Dumux
{
// \{

/*!
 * \brief Labels for pores and throats
 */
struct Labels
{
    static const int interior = -1; //!< Label for pores/throats not on a boundary
    static const int noflow = 0; //!< Label for pores/throats with a now flow BC
    static const int dirichlet = 1; //!< Label for pores/throats with Dirichlet BC
    static const int inlet = 2; //!< Label for pores/throats on an inlet
    static const int outlet = 3; //!< Label for pores/throats on an outlet
    static const int source = 4; //!< Label for pores/throats with a sink/source term

    static const int coupled1 = 11; //!< Label for pores/throats coupled to another model
    static const int coupled2 = 12; //!< Label for pores/throats coupled to another model
};

// \}
} // end namespace

#endif
