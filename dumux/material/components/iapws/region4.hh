// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup IAPWS
 * \brief Implements the equations for region 4 of the IAPWS '97 formulation.
 *
 * See:
 *
 * IAPWS: "Revised Release on the IAPWS Industrial Formulation
 * 1997 for the Thermodynamic Properties of Water and Steam",
 * http://www.iapws.org/relguide/IF97-Rev.pdf
 */
#ifndef DUMUX_IAPWS_REGION4_HH
#define DUMUX_IAPWS_REGION4_HH

#include <cmath>
#include <iostream>

#include <dune/common/math.hh>

namespace Dumux::IAPWS {

/*!
 * \ingroup IAPWS
 * \brief Implements the equations for region 4 of the IAPWS '97 formulation.
 * \tparam Scalar The type used for scalar values
 *
 * See:
 *
 * IAPWS: "Revised Release on the IAPWS Industrial Formulation
 * 1997 for the Thermodynamic Properties of Water and Steam",
 * http://www.iapws.org/relguide/IF97-Rev.pdf
 */
template <class Scalar>
class Region4
{
public:
    /*!
     * \brief Returns the saturation pressure in \f$\mathrm{[Pa]}\f$ of pure water at a given
     *        temperature.
     *
     *\param temperature temperature of component in \f$\mathrm{[K]}\f$
     *
     * The saturation pressure is often also called vapor pressure. This formulation is valid
     * for a temperature range between 273.15 K and 647.096 K (critical temperature).
     */
    static Scalar saturationPressure(Scalar temperature)
    {
        constexpr Scalar n[10] = {
            0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
            0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
            -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
            0.65017534844798e3
        };

        const Scalar sigma = temperature + n[8]/(temperature - n[9]);

        const Scalar A = (sigma + n[0])*sigma + n[1];
        const Scalar B = (n[2]*sigma + n[3])*sigma + n[4];
        const Scalar C = (n[5]*sigma + n[6])*sigma + n[7];

        using std::sqrt;
        Scalar tmp = 2*C/(sqrt(B*B - 4*A*C) - B);
        tmp *= tmp;
        tmp *= tmp;

        return 1e6*tmp;
    }

    /*!
     * \brief Returns the saturation temperature in \f$\mathrm{[K]}\f$ of pure water at a given
     *        pressure.
     *
     *\param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * The saturation pressure is often also called vapor pressure.
     */
    static Scalar vaporTemperature(Scalar pressure)
    {
        constexpr Scalar n[10] = {
            0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
            0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
            -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
            0.65017534844798e3
        };

        using std::pow;
        using Dune::power;
        const Scalar beta = pow((pressure/1e6 /*from Pa to MPa*/), (1./4.));
        const Scalar beta2 = power(beta, 2);
        const Scalar E = beta2 + n[2] * beta + n[5];
        const Scalar F = n[0]*beta2 + n[3]*beta + n[6];
        const Scalar G = n[1]*beta2 + n[4]*beta + n[7];

        using std::sqrt;
        const Scalar D = ( 2.*G)/(-F -sqrt(power(F,2) - 4.*E*G));
        return (n[9] + D - sqrt(power(n[9]+D , 2) - 4.* (n[8] + n[9]*D)) ) * 0.5;
    }
};

} // end namespace Dumux::IAPWS

#endif
