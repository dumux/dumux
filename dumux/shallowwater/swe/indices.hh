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
 * \ingroup SweModel
 * \copydoc Dumux::SweIndices
 */
#ifndef DUMUX_SWE_INDICES_HH
#define DUMUX_SWE_INDICES_HH

#include <dumux/common/properties.hh>

namespace Dumux
{
// \{
/*!
 * \ingroup SweModel
 * \brief The common indices for the shallow water equations model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
struct SweIndices
{

    static constexpr int dimXIdx = 0; //!< Index of the x-component of a vector of size dim
    static constexpr int dimYIdx = 1; //!< Index of the y-component of a vector of size dim

    static constexpr int massBalanceIdx = 0; //!< Index of the mass balance equation
    static constexpr int momentumXBalanceIdx = 1; //!< Index of the x momentum balance equation
    static constexpr int momentumYBalanceIdx = 2; //!< Index of the y momentum balance equation

    static constexpr auto dim = 2; //!< We have allways two-dimensional SWEs
    static constexpr auto numEq = 3;  //!< SWEs have allways 3 balance equations for mass and momentum

    static constexpr int waterdepthIdx = massBalanceIdx; //!< Index of the velocity in a solution vector
    static constexpr int velocityXIdx = momentumXBalanceIdx; //!< Index of the velocity in a solution vector
    static constexpr int velocityYIdx = momentumYBalanceIdx; //!< Index of the velocity in a solution vector

};

// \}
} // end namespace

#endif
