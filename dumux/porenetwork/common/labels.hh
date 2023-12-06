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
 * \ingroup PoreNetworkModels
 * \brief Defines labels for pores and throats.
 */
#ifndef DUMUX_PNM_LABELS_HH
#define DUMUX_PNM_LABELS_HH

namespace Dumux::PoreNetwork {
// \{

/*!
 * \ingroup PoreNetworkModels
 * \brief Labels for pores and throats.
 */
struct Labels
{
    // TODO revise this concept
    static constexpr int interior = -1; //!< Label for pores/throats not on a boundary
    static constexpr int noflow = 0; //!< Label for pores/throats with a now flow BC
    static constexpr int dirichlet = 1; //!< Label for pores/throats with Dirichlet BC
    static constexpr int inlet = 2; //!< Label for pores/throats on an inlet
    static constexpr int outlet = 3; //!< Label for pores/throats on an outlet
    static constexpr int source = 4; //!< Label for pores/throats with a sink/source term
    static constexpr int internalDirichlet = 10; //!< Label for internal pores/throats with Dirichlet BC
};

// \}
} // end Dumux::PoreNetwork

#endif
