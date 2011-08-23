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

extern "C" {
#include <quadmath.h>
}

typedef __float128 quad;

namespace std
{

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
