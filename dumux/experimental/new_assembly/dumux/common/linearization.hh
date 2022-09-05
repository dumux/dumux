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
 * \copydoc Dumux::Linearization
 */
#ifndef DUMUX_COMMON_LINEARIZATION_HH
#define DUMUX_COMMON_LINEARIZATION_HH

#include <type_traits>
#include <concepts>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief Class to represent a linearized function, storing
 *        references to the derivative and function value.
 */
template<typename D, typename V>
class Linearization
{
    template<typename Arg, typename Member>
    static constexpr bool isValidConstructorArg = std::is_lvalue_reference_v<Arg> &&
                                                  std::convertible_to<Arg, Member&>;

public:
    using Derivative = D;
    using Value = V;

    /*!
    * \brief Constructor from derivative and function value
    * \note This is templated such that we can check at compile-time that
    *       l-value references are passed in. If we made this simply const
    *       refs, then one may pass in temporaries and cause segfaults.
    */
    template<typename _D, typename _R> requires(
        isValidConstructorArg<_D, const Derivative> and
        isValidConstructorArg<_R, const V>)
    Linearization(_D&& derivative, _R&& value)
    : derivative_(derivative)
    , value_(value)
    {}

    const Derivative& derivative() const { return derivative_; }
    const Value& value() const { return value_; }

private:
    const Derivative& derivative_;
    const Value& value_;
};

template<typename D, typename R>
Linearization(D&&, R&&) -> Linearization<std::decay_t<D>, std::decay_t<R>>;

} // end namespace Dumux

#endif
