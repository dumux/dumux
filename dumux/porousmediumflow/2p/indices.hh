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
 * \ingroup TwoPModel
 * \brief Defines the indices required for the two-phase fully implicit model.
 */
#ifndef DUMUX_BOX_2P_INDICES_HH
#define DUMUX_BOX_2P_INDICES_HH

#include <dumux/common/properties.hh>
#include "formulation.hh"

namespace Dumux
{
// \{

/*!
 * \ingroup TwoPModel
 * \brief Defines the indices required for the two-phase fully implicit model.
 *
 * \tparam FluidSystem The fluid system class
 * \tparam PVOffset The first index in a primary variable vector.
 */
struct TwoPCommonIndices
{
    // Primary variable indices
    static const int pressureIdx = 0; //!< index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int saturationIdx = 1; //!< index of the saturation of the non-wetting/wetting phase

    // indices of the equations
    static const int conti0EqIdx = 0; //!< index of the first continuity equation
    static const int contiWEqIdx = 0; //!< index of the continuity equation of the wetting phase
    static const int contiNEqIdx = 1; //!< index of the continuity equation of the non-wetting phase
};


/*!
 * \ingroup TwoPModel
 * \brief The indices for the \f$p_w-S_n\f$ formulation of the
 *        isothermal two-phase model.
 *
 * \tparam FluidSystem The fluid system class
 * \tparam formulation The formulation, either pwsn or pnsw
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <TwoPFormulation formulation = TwoPFormulation::pwsn>
struct TwoPIndices : public TwoPCommonIndices
{
    // indices of the primary variables
    static constexpr int pwIdx = 0; //!< index of the wetting phase pressure
    static constexpr int snIdx = 1; //!< index of the nonwetting phase saturation
};

/*!
 * \ingroup TwoPModel
 * \brief The indices for the \f$p_n-S_w\f$ formulation of the
 *        isothermal two-phase model.
 *
 * \tparam FluidSystem The fluid system class
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <>
struct TwoPIndices<TwoPFormulation::pnsw>
: public TwoPCommonIndices
{
    // indices of the primary variables
    static constexpr int pnIdx = 0; //!< index of the nonwetting phase pressure
    static constexpr int swIdx = 1; //!< index of the wetting phase saturation
};

// \}
} // namespace Dumux


#endif
