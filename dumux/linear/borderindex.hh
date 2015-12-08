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
 *
 * \brief A single index intersecting with the process boundary.
 */
#ifndef DUMUX_BORDER_INDEX_HH
#define DUMUX_BORDER_INDEX_HH

#warning This file is deprecated and will be removed after Dumux 2.9

namespace Dumux {

/*!
 * \brief A single index intersecting with the process boundary.
 */
struct BorderIndex
{
    //! Index of the entity for the local process
    int localIdx;

    //! Index of the entity for the peer process
    int peerIdx;

    //! Rank of the peer process
    int peerRank;

    //! Distance to the process border for the peer (in hops)
    int borderDistance;

    //! True if and only if the entity which corrosponds to the index
    //! is in the interior of more than one process (i.e. the entity
    //! is completely on the border)
    bool isShared;
};

}

#endif
