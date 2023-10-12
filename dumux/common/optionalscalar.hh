// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief A wrapper that can either contain a valid Scalar or NaN
 */
#ifndef DUMUX_COMMON_OPTIONAL_SCALAR_HH
#define DUMUX_COMMON_OPTIONAL_SCALAR_HH

#include <limits>
#include <cmath>

namespace Dumux {

/*!
 * \ingroup Core
 * \brief A type for an optional scalar (contains either a valid number or NaN)
 * \tparam T Type of the underlying floating point number type
 */
template<class T>
struct OptionalScalar
{
    static_assert(std::numeric_limits<T>::has_quiet_NaN, "T has to be able to represent a quiet NaN!");

    OptionalScalar() = default;

    OptionalScalar(T value)
    : value_(value)
    {}

    T value() const
    { return value_; }

    explicit operator bool() const
    {
        using std::isnan;
        return !isnan(value_);
    }
private:
    T value_ = std::numeric_limits<T>::quiet_NaN();
};

} // namespace Dumux

#endif
