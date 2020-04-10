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
 * \ingroup Elastic
 * \brief Defines the indices for the elastic model
 */
#ifndef DUMUX_ELASTIC_INDICES_HH
#define DUMUX_ELASTIC_INDICES_HH

namespace Dumux {

/*!
 * \ingroup Elastic
 * \brief The indices for the linear elasticity model.
 */
struct ElasticIndices
{
    // returns the equation index for a given space direction
    static constexpr int momentum(int dirIdx) { return dirIdx; };

    // returns the primary variable index for a given space direction
    static constexpr int u(int dirIdx) { return dirIdx; };

    // Equation indices
    static constexpr int momentumXEqIdx = 0;
    static constexpr int momentumYEqIdx = 1;
    static constexpr int momentumZEqIdx = 2;

    // primary variable indices
    static constexpr int uxIdx = 0;
    static constexpr int uyIdx = 1;
    static constexpr int uzIdx = 2;
};

} // end namespace Dumux

#endif
