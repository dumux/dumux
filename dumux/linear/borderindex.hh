/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
 * \brief Represents a single index intersecting with the process
 *        boundary.
 */
#ifndef DUMUX_BORDER_INDEX_HH
#define DUMUX_BORDER_INDEX_HH

namespace Dumux {

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
};

};

#endif
