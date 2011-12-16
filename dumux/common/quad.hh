// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * \brief This file provides the infrastructure to use quad-precision
 *        floating point values in the numerical models.
 */
#ifndef DUMUX_QUAD_HH
#define DUMUX_QUAD_HH

#if ! HAVE_QUAD
#error "Quad precision floating point values must be enabled to use this file"
#endif

#include <iostream>
#include <cmath>
#include <limits>

extern "C" {
#include <quadmath.h>
}

typedef __float128 quad;

namespace std
{

// provide the numeric limits for the quad precision type
template <>
class numeric_limits<quad>
{
public:
    static constexpr bool is_specialized = true;

    static constexpr quad min() throw()
    { return FLT128_MIN; }
    static constexpr quad max() throw()
    { return FLT128_MAX; }

    // number of bits in mantissa
    static constexpr int  digits = FLT128_MANT_DIG;
    // number of decimal digits
    static constexpr int  digits10 = FLT128_DIG;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = false;
    static constexpr bool is_exact = false;
    static constexpr int radix = 0;
    static constexpr quad epsilon() throw()
    { return FLT128_EPSILON; }
    static constexpr quad round_error() throw()
    { return 0.5; }

    static constexpr int min_exponent = FLT128_MIN_EXP;
    static constexpr int min_exponent10 = FLT128_MIN_10_EXP;
    static constexpr int max_exponent = FLT128_MAX_EXP;
    static constexpr int max_exponent10 = FLT128_MAX_10_EXP;

    static constexpr bool has_infinity = true;
    static constexpr bool has_quiet_NaN = true;
    static constexpr bool has_signaling_NaN = true;
    static constexpr float_denorm_style has_denorm = denorm_present;
    static constexpr bool has_denorm_loss = false;
    static constexpr quad infinity() throw()
    { return __builtin_huge_valq(); };
    static constexpr quad quiet_NaN() throw()
    { return __builtin_nan(""); }
    static constexpr quad signaling_NaN() throw()
    { return __builtin_nans(""); }
    static constexpr quad denorm_min() throw()
    { return FLT128_DENORM_MIN; }

    static constexpr bool is_iec559 = true;
    static constexpr bool is_bounded = true;
    static constexpr bool is_modulo = false;

    static constexpr bool traps = std::numeric_limits<double>::traps;
    static constexpr bool tinyness_before = std::numeric_limits<double>::tinyness_before;
    static constexpr float_round_style round_style = round_to_nearest;
};

inline std::ostream& operator<<(std::ostream& os, const quad &val)
{
    return (os << double(val));
}

inline std::istream& operator>>(std::istream& is, quad &val)
{
    double tmp;
    std::istream &ret = (is >> tmp);
    val = tmp;
    return ret;
}

inline quad abs(quad val)
{ return (val < 0)?-val:val; }

inline quad max(quad a, quad b)
{ return (a>b)?a:b; };

inline quad min(quad a, quad b)
{ return (a<b)?a:b; };

inline quad sqrt(quad val)
{ return sqrtq(val); };

template <class ExpType>
inline quad pow(quad base, ExpType exp)
{ return powq(base, exp); };

inline quad exp(quad val)
{ return expq(val); };

inline bool isfinite(quad val)
{ return !isnanq(val) && !isinfq(val); };

inline bool isnan(quad val)
{ return isnanq(val); };

inline bool isinf(quad val)
{ return isinfq(val); };

} // namespace std


#endif // DUMUX_QUAD_HH
