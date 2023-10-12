// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief The cubic hermite spline basis
 */
#ifndef DUMUX_COMMON_CUBIC_SPLINE_HERMITE_BASIS_HH
#define DUMUX_COMMON_CUBIC_SPLINE_HERMITE_BASIS_HH

namespace Dumux {

/*!
 * \ingroup Core
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
