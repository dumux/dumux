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
 * \ingroup Common
 * \brief This file provides the infrastructure to use quad-precision
 *        floating point values in the numerical models.
 */
#if !defined DUMUX_QUAD_HH && HAVE_QUAD
#define DUMUX_QUAD_HH

#include <iostream>
#include <cmath>
#include <limits>

extern "C" {
#include <quadmath.h>
}

#include <dune/common/typetraits.hh>

namespace Dumux {

using quad = __float128;

} // namespace Dumux

// Dune's desired way of enabling their algebraic operations for
// extended-precision data types, see dune/common/typetraits.hh.
namespace Dune {

template <>
struct IsNumber<Dumux::quad> : std::true_type {};

} // namespace Dune

namespace std {

// provide the numeric limits for the quad precision type
template <>
class numeric_limits<Dumux::quad>
{
public:
    static constexpr bool is_specialized = true;

    static constexpr Dumux::quad min() noexcept
    { return FLT128_MIN; }
    static constexpr Dumux::quad max() noexcept
    { return FLT128_MAX; }
    static constexpr Dumux::quad lowest() noexcept
    { return -FLT128_MAX; }

    // number of bits in mantissa
    static constexpr int digits = FLT128_MANT_DIG;
    // number of decimal digits
    static constexpr int digits10 = FLT128_DIG;
    static constexpr int max_digits10 = FLT128_MANT_DIG;

    static constexpr bool is_signed = true;
    static constexpr bool is_integer = false;
    static constexpr bool is_exact = false;
    static constexpr int radix = 0;
    static constexpr Dumux::quad epsilon() noexcept
    { return FLT128_EPSILON; }
    static constexpr Dumux::quad round_error() noexcept
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
    static constexpr Dumux::quad infinity() noexcept
    { return __builtin_huge_valq(); };
    static constexpr Dumux::quad quiet_NaN() noexcept
    { return __builtin_nan(""); }
    static constexpr Dumux::quad signaling_NaN() noexcept
    { return __builtin_nans(""); }
    static constexpr Dumux::quad denorm_min() noexcept
    { return FLT128_DENORM_MIN; }

    static constexpr bool is_iec559 = true;
    static constexpr bool is_bounded = true;
    static constexpr bool is_modulo = false;

    static constexpr bool traps = std::numeric_limits<double>::traps;
    static constexpr bool tinyness_before = std::numeric_limits<double>::tinyness_before;
    static constexpr float_round_style round_style = round_to_nearest;
};

} // namespace std

// Putting the following definitions in the namespace std yields undefined
// behavior. Unfortunately, ADL doesn't work for aliases (Quad) because the underlying type
// is used for the lookup. Putting these functions into any other namespace yields compiler errors.
namespace std {

inline ostream& operator<<(ostream& os, const Dumux::quad &val)
{
    return (os << double(val));
}

inline istream& operator>>(istream& is, Dumux::quad &val)
{
    double tmp;
    istream &ret = (is >> tmp);
    val = tmp;
    return ret;
}

inline Dumux::quad abs(Dumux::quad val)
{ return (val < 0)?-val:val; }

inline Dumux::quad floor(Dumux::quad val)
{ return floorq(val); }

inline Dumux::quad ceil(Dumux::quad val)
{ return ceilq(val); }

inline Dumux::quad max(Dumux::quad a, Dumux::quad b)
{ return (a>b)?a:b; }

inline Dumux::quad min(Dumux::quad a, Dumux::quad b)
{ return (a<b)?a:b; }

inline Dumux::quad sqrt(Dumux::quad val)
{ return sqrtq(val); }

template <class ExpType>
inline Dumux::quad pow(Dumux::quad base, ExpType exp)
{ return powq(base, exp); }

inline Dumux::quad exp(Dumux::quad val)
{ return expq(val); }

inline Dumux::quad log(Dumux::quad val)
{ return logq(val); }

inline Dumux::quad sin(Dumux::quad val)
{ return sinq(val); }

inline Dumux::quad cos(Dumux::quad val)
{ return cosq(val); }

inline Dumux::quad tan(Dumux::quad val)
{ return tanq(val); }

inline Dumux::quad atan(Dumux::quad val)
{ return atanq(val); }

inline Dumux::quad atan2(Dumux::quad a, Dumux::quad b)
{ return atan2q(a, b); }

inline bool signbit(Dumux::quad val)
{ return signbitq(val); }

inline bool isfinite(Dumux::quad val)
{ return !isnanq(val) && !isinfq(val); }

inline bool isnan(Dumux::quad val)
{ return isnanq(val); }

inline bool isinf(Dumux::quad val)
{ return isinfq(val); }

} // namespace std

#endif // DUMUX_QUAD_HH
