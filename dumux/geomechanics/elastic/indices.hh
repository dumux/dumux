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
 * \brief Defines the primary variable and equation indices used by
 *        the linear elasticity model
 */
#ifndef DUMUX_ELASTIC_INDICES_HH
#define DUMUX_ELASTIC_INDICES_HH

namespace Dumux
{
// \{

/*!
 * \ingroup ElasticBoxModel
 * \ingroup ImplicitIndices
 * \brief The indices for the linear elasticity model.
 */
template <int PVOffset = 0>
struct ElasticIndices
{
    // returns the equation index for a given space direction
    static int momentum(int dirIdx)
    {
        return PVOffset + dirIdx;
    };

    // returns the primary variable index for a given space direction
    static int u(int dirIdx)
    {
        return PVOffset + dirIdx;
    };

    // Equation indices
    static const int momentumXEqIdx = PVOffset + 0;
    static const int momentumYEqIdx = PVOffset + 1;
    static const int momentumZEqIdx = PVOffset + 2;

    // primary variable indices
    static const int uxIdx = PVOffset + 0;
    static const int uyIdx = PVOffset + 1;
    static const int uzIdx = PVOffset + 2;
};

}// namespace Dumux

#endif
