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
 * \brief Binary coefficients for water and methane.
 */
#ifndef DUMUX_BINARY_COEFF_H2O_CH4_HH
#define DUMUX_BINARY_COEFF_H2O_CH4_HH

#include <dumux/material/binarycoefficients/henryiapws.hh>
#include <dumux/material/binarycoefficients/fullermethod.hh>

#include <dumux/material/components/ch4.hh>
#include <dumux/material/components/h2o.hh>

namespace Dumux {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and methane.
 */
class H2O_CH4
{
public:
    /*!
     * \brief Henry coefficient \f$[N/m^2]\f$  for molecular methane in liquid water.
     *
     * See:
     *
     * IAPWS: "Guideline on the Henry's Constant and Vapor-Liquid
     * Distribution Constant for Gases in H2O and D2O at High
     * Temperatures"
     * http://www.iapws.org/relguide/HenGuide.pdf
     */
    template <class Scalar>
    static Scalar henry(Scalar temperature)
    {
        const Scalar E = 2215.6977;
        const Scalar F = -0.1089;
        const Scalar G = -6.6240;
        const Scalar H = 4.6789;

        return henryIAPWS(E, F, G, H, temperature);
    }

    /*!
     * \brief Binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular water in methane.
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        //         DUNE_THROW(Dune::NotImplemented, "diffusion coefficient for gasous water and methane");

        typedef Dumux::Components::H2O<Scalar> H2O;
        typedef Dumux::Components::CH4<Scalar> CH4;

        // atomic diffusion volumes
        // Vch4 = sum(ni*Vi) = 15.9 + 4*2.31 = 25.14 (Tang et al., 2015)--> method, (Poling et al., 2001, p.11.11)--> values
        const Scalar SigmaNu[2] = { 13.1 /* H2O */,  25.14 /* CH4 */ };
        // molar masses [g/mol]
        const Scalar M[2] = { H2O::molarMass()*Scalar(1e3), CH4::molarMass()*Scalar(1e3) };

        return fullerMethod(M, SigmaNu, temperature, pressure);
    }

    /*!
     * \brief Diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular methane in liquid water.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     *
     * The empirical equations for estimating the diffusion coefficient in
     * infinite solution which are presented in Reid, 1987 \cite reid1987 all show a
     * linear dependency on temperature. We thus simply scale the
     * experimentally obtained diffusion coefficient of Ferrell and
     * Himmelblau by the temperature.<br>
     * This function use an interpolation of the data by \cite witherspoon1965
     * http://dx.doi.org/10.1021/j100895a017
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        return 2.93856e-11 * temperature - 6.89402e-09;
    }
};

} // end namespace BinaryCoeff
} // end namespace Dumux

#endif
