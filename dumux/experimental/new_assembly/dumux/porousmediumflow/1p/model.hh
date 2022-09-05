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
 * \ingroup PorousMediumFlow
 * \brief Model for single-phase flow in porous-media.
 */
#ifndef DUMUX_POROUS_MEDIUM_FLOW_ONEP_MODEL_HH
#define DUMUX_POROUS_MEDIUM_FLOW_ONEP_MODEL_HH

#include <utility>
#include <concepts>
#include <type_traits>

namespace Dumux {

template<typename Scalar,
         typename FluidSystem>
class OnePModel
{
public:
    static constexpr int numEq = 1;

    using VolumeVariables = SomeVolVars<FluidSystem, ...>;

    struct Indices
    {
        static constexpr int pressureIdx = 0;
        static constexpr int conti0EqIdx = 0;
    };

    auto storageOperator() const
    { ... }

    auto fluxOperator() const
    {
        // TODO: Some flux thing that sums up all fluxes?
        // But also allow extracting single ones? (see advectiveFluxoperator())
        // Also, the flux thing should be as efficient as possible, e.g. by caching
        // KgradP once and then only apply different upwind terms...
        ...
    }

    auto advectiveFluxOperator() const { ... }
};

} // namespace Dumux

#endif
