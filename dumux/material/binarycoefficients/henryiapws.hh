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
 * \ingroup Binarycoefficients
 * \brief The IAPWS formulation of Henry coefficients in water.
 */
#ifndef DUMUX_HENRY_IAPWS_HH
#define DUMUX_HENRY_IAPWS_HH

#include <dumux/material/components/h2o.hh>

namespace Dumux {
/*!
 * \ingroup Binarycoefficients
 * \brief The Henry constants in liquid water using the IAPWS 2004
 *        formulation.
 * \param E Correlation parameter
 * \param F Correlation parameter
 * \param G Correlation parameter
 * \param H Correlation parameter
 * \param temperature the temperature \f$\mathrm{[K]}\f$
 *
 * This function calculates \f$\mathrm{K_D}\f$, see:
 *
 * IAPWS: "Guideline on the Henry's Constant and Vapor-Liquid
 * Distribution Constant for Gases in H2O and D2O at High
 * Temperatures"
 * Equation (5) \cite watanabe2004 <BR>
 *
 * Range of validity: T = {278.12 ; 636.46}
 * approximations beyond this range are increasingly incorrect.
 * However, close to the critical the values are more, again.
 */
template <class Scalar>
inline Scalar henryIAPWS(Scalar E,
                         Scalar F,
                         Scalar G,
                         Scalar H,
                         Scalar temperature)
{
    using H2O = Dumux::Components::H2O<Scalar>;

    // regularizing temperature helps for stability.
    // Results are unphysical!
    if (temperature > H2O::criticalTemperature())
        temperature = H2O::criticalTemperature();

    Scalar Tr = temperature/H2O::criticalTemperature();
    Scalar tau = 1 - Tr;

    static const Scalar c[6] = {
        1.99274064, 1.09965342, -0.510839303,
        -1.75493479,-45.5170352, -6.7469445e5
    };
    static const Scalar d[6] = {
        1/3.0, 2/3.0, 5/3.0,
        16/3.0, 43/3.0, 110/3.0
    };
    static const Scalar q = -0.023767;

    Scalar f = 0;
    using std::pow;
    for (int i = 0; i < 6; ++i) {
        f += c[i]*pow(tau, d[i]);
    }

    Scalar exponent =
        q*F +
        E/temperature*f +
        (F +
         G*pow(tau, 2.0/3) +
         H*tau)*
        exp((H2O::tripleTemperature() - temperature)/100);

    using std::exp;
    return exp(exponent)*H2O::vaporPressure(temperature);
}
}

#endif
