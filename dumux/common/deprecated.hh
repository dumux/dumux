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
 * \ingroup Common
 * \brief Helpers for deprecation
 */

#ifndef DUMUX_COMMON_DEPRECATED_HH
#define DUMUX_COMMON_DEPRECATED_HH

#include <utility>

#include <dune/common/ftraits.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/std/type_traits.hh>
#include <dumux/common/numeqvector.hh>

namespace Dumux {

#ifndef DOXYGEN // hide from doxygen
// Helper classes/functions for deprecation
// Each implementation has to state after which release
// it will be removed. Implementations in the Deprecated
// namespace will be removed without
// deprecation after their usage in the code expired,
// so most likely you don't want to use this in your code
namespace Deprecated {

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif // __clang__

#ifdef __clang__
#pragma clang diagnostic pop
#endif  // __clang__

template<class CO2Impl>
    struct BrineCO2Helper{
        template<class CO2Arg>
        using TabulatedDensityDetector = decltype(std::declval<CO2Arg>().tabulatedDensity);
        static constexpr bool rawCO2Table = Dune::Std::is_detected<TabulatedDensityDetector,
                                                                   CO2Impl >::value;

        template< typename T>
        [[deprecated("Passing just CO2Tables to define a BrineCO2 fluidsystem/binarycoefficient is deprecated. Use Components::CO2<Scalar, CO2Tables> as template parameter instead.")]]
        static constexpr void DefiningBrineCO2WithCO2Table() {}

        static constexpr bool isRawTable()
        {
            if constexpr (rawCO2Table)
                DefiningBrineCO2WithCO2Table<CO2Impl>();
            return rawCO2Table;
        }
    };

} // end namespace Deprecated
#endif

} // end namespace Dune
#endif
