// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Deflection of an annular Kirchhoff-Love plate, clamped inside and free outside
 */
#ifndef DUMUX_TEST_SOLIDMECHANICS_PLATE_ANNULUS_PLATE_HH
#define DUMUX_TEST_SOLIDMECHANICS_PLATE_ANNULUS_PLATE_HH

#include <cmath>
#include <dune/common/fvector.hh>

namespace Dumux {

/*!
 * \brief Coefficients of the axisymmetric deflection of an annular plate that is
 *        clamped at \f$ r=a \f$ and free at \f$ r=b \f$ under a uniform load
 *
 * The general axisymmetric solution of \f$ -D\Delta^2 w = q \f$ is
 * \f$ w = -qr^4/(64D) + C_1 + C_2 r^2 + C_3\ln r + C_4 r^2\ln r \f$. The free edge
 * conditions \f$ M_{rr}(b) = 0 \f$ and \f$ Q_r(b) = 0 \f$ fix \f$ C_4 \f$ and one
 * relation between \f$ C_2 \f$ and \f$ C_3 \f$, the clamped conditions
 * \f$ w(a) = w'(a) = 0 \f$ the remaining two.
 */
template<class Scalar>
Dune::FieldVector<Scalar, 4>
annulusPlateCoefficients(Scalar a, Scalar b, Scalar q, Scalar D, Scalar nu)
{
    using std::log;
    const auto c4 = q*b*b/(8.0*D);
    const auto rhsMoment = q*b*b*(3.0 + nu)/(16.0*D) - c4*(2.0*(1.0 + nu)*log(b) + 3.0 + nu);
    const auto rhsSlope = q*a*a*a/(16.0*D) - c4*a*(2.0*log(a) + 1.0);
    const auto c2 = (rhsMoment + (1.0 - nu)*a*rhsSlope/(b*b))
                    /(2.0*(1.0 + nu) + 2.0*a*a*(1.0 - nu)/(b*b));
    const auto c3 = a*rhsSlope - 2.0*c2*a*a;
    const auto c1 = q*a*a*a*a/(64.0*D) - c2*a*a - c3*log(a) - c4*a*a*log(a);
    return {c1, c2, c3, c4};
}

//! Axisymmetric deflection at radius r with coefficients C1, C2, C3, C4
template<class Scalar, class Coefficients>
Scalar annulusPlateDeflection(Scalar r, Scalar q, Scalar D, const Coefficients& c)
{
    using std::log;
    return -q*r*r*r*r/(64.0*D) + c[0] + c[1]*r*r + c[2]*log(r) + c[3]*r*r*log(r);
}

} // end namespace Dumux

#endif
