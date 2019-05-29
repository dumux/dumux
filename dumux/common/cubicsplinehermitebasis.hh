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
 * \brief The cubic hermite spline basis
 */
#ifndef DUMUX_COMMON_CUBIC_SPLINE_HERMITE_BASIS_HH
#define DUMUX_COMMON_CUBIC_SPLINE_HERMITE_BASIS_HH

namespace Dumux {

/*!
 * \ingroup Common
 * \brief The cubic spline hermite basis
 * \note see https://en.wikipedia.org/wiki/Cubic_Hermite_spline
 */
template<class Scalar = double>
struct CubicSplineHermiteBasis
{
    static constexpr Scalar h00(const Scalar t)
    { return t*(2.0*t*t - 3.0*t) + 1.0; }

    static constexpr Scalar h10(const Scalar t)
    { return t*(t*t - 2.0*t + 1.0); }

    static constexpr Scalar h01(const Scalar t)
    { return t*t*(3.0 - 2.0*t); }

    static constexpr Scalar h11(const Scalar t)
    { return t*t*(t - 1.0); }

    static constexpr Scalar dh00(const Scalar t)
    { return 6.0*t*(t - 1.0); }

    static constexpr Scalar dh10(const Scalar t)
    { return t*(3.0*t - 4.0) + 1.0; }

    static constexpr Scalar dh01(const Scalar t)
    { return 6.0*t*(1.0 - t); }

    static constexpr Scalar dh11(const Scalar t)
    { return t*(3.0*t - 2.0); }
};

} // end namespace Dumux

#endif
