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
 * \brief A wrapper that can either contain an object of T or be empty.
 */
#ifndef DUMUX_COMMON_OPTIONAL_HH
#define DUMUX_COMMON_OPTIONAL_HH

#include <utility>
#include <limits>
#include <cmath>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief A type for an optional floating point number (if NaN can be spared as a valid result)
 * \tparam T Type of wrapped floating point number
 */
template<class T>
struct OptionalFloatingPointNumber
{
    static_assert(std::numeric_limits<T>::has_quiet_NaN, "T has to be able to represent a quiet NaN!");

    OptionalFloatingPointNumber() = default;

    OptionalFloatingPointNumber(T value) : value_(value) {}

    T value() const
    { return value_; }

    operator bool() const
    {
        using std::isnan;
        return !isnan(value_);
    }
private:
    T value_ = std::numeric_limits<T>::quiet_NaN();
};

} // namespace Dumux

#endif // DUMUX_COMMON_OPTIONAL_HH
