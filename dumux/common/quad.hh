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

using Quad = __float128;

} // namespace Dumux

// Dune's desired way of enabling their algebraic operations for
// extended-precision data types, see dune/common/typetraits.hh.
namespace Dune {

template <>
struct IsNumber<Dumux::Quad> : std::true_type {};

} // namespace Dune

namespace std {

// provide the numeric limits for the quad precision type
template <>
class numeric_limits<Dumux::Quad>
{
public:
    static constexpr bool is_specialized = true;

    static constexpr Dumux::Quad min() noexcept
    { return FLT128_MIN; }
    static constexpr Dumux::Quad max() noexcept
    { return FLT128_MAX; }
    static constexpr Dumux::Quad lowest() noexcept
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
    static constexpr Dumux::Quad epsilon() noexcept
    { return FLT128_EPSILON; }
    static constexpr Dumux::Quad round_error() noexcept
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
    static constexpr Dumux::Quad infinity() noexcept
    { return __builtin_huge_valq(); };
    static constexpr Dumux::Quad quiet_NaN() noexcept
    { return __builtin_nan(""); }
    static constexpr Dumux::Quad signaling_NaN() noexcept
    { return __builtin_nans(""); }
    static constexpr Dumux::Quad denorm_min() noexcept
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

inline ostream& operator<<(ostream& os, const Dumux::Quad &val)
{
    return (os << double(val));
}

inline istream& operator>>(istream& is, Dumux::Quad &val)
{
    double tmp;
    istream &ret = (is >> tmp);
    val = tmp;
    return ret;
}

inline Dumux::Quad abs(Dumux::Quad val)
{ return (val < 0)?-val:val; }

inline Dumux::Quad floor(Dumux::Quad val)
{ return floorq(val); }

inline Dumux::Quad ceil(Dumux::Quad val)
{ return ceilq(val); }

inline Dumux::Quad max(Dumux::Quad a, Dumux::Quad b)
{ return (a>b)?a:b; }

inline Dumux::Quad min(Dumux::Quad a, Dumux::Quad b)
{ return (a<b)?a:b; }

inline Dumux::Quad sqrt(Dumux::Quad val)
{ return sqrtq(val); }

template <class ExpType>
inline Dumux::Quad pow(Dumux::Quad base, ExpType exp)
{ return powq(base, exp); }

inline Dumux::Quad exp(Dumux::Quad val)
{ return expq(val); }

inline Dumux::Quad log(Dumux::Quad val)
{ return logq(val); }

inline Dumux::Quad sin(Dumux::Quad val)
{ return sinq(val); }

inline Dumux::Quad cos(Dumux::Quad val)
{ return cosq(val); }

inline Dumux::Quad tan(Dumux::Quad val)
{ return tanq(val); }

inline Dumux::Quad atan(Dumux::Quad val)
{ return atanq(val); }

inline Dumux::Quad atan2(Dumux::Quad a, Dumux::Quad b)
{ return atan2q(a, b); }

inline bool signbit(Dumux::Quad val)
{ return signbitq(val); }

inline bool isfinite(Dumux::Quad val)
{ return !isnanq(val) && !isinfq(val); }

inline bool isnan(Dumux::Quad val)
{ return isnanq(val); }

inline bool isinf(Dumux::Quad val)
{ return isinfq(val); }

} // namespace std

#endif // DUMUX_QUAD_HH
