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
 * \brief Binary coefficients for air and xylene.
 */
#ifndef DUMUX_BINARY_COEFF_AIR_XYLENE_HH
#define DUMUX_BINARY_COEFF_AIR_XYLENE_HH

#include <algorithm>

#include <dune/common/math.hh>

#include <dumux/material/components/air.hh>
#include <dumux/material/components/xylene.hh>

namespace Dumux::BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for air and xylene.
 */
class Air_Xylene
{
public:
    /*!
     * \brief Henry coefficient \f$\mathrm{[Pa]}\f$  for mesitylene in air.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     */
    template <class Scalar>
    static Scalar henry(Scalar temperature)
    { DUNE_THROW(Dune::NotImplemented,
                 "Henry coefficient of air in xylene");
    }

    /*!
     * \brief Binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for air and xylene.
     * method according to Wilke and Lee
     * see W.J. Lyman, W.F. Reehl, D.H. Rosenblatt (1990) \cite lyman1990 <BR>
     * \param temperature temperature in \f$\mathrm{[K]}\f$
     * \param pressure pressure in \f$\mathrm{[Pa]}\f$
     *
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        using Air = Dumux::Components::Air<Scalar>;
        using Xylene = Dumux::Components::Xylene<Scalar>;

        using std::clamp;
        temperature = clamp(temperature, 1e-9, 500.0); // regularization
        pressure = clamp(pressure, 0.0, 1e8); // regularization

        using std::pow;
        using Dune::power;
        using std::sqrt;
        using std::exp;
        const Scalar M_x = 1e3*Xylene::molarMass(); // [g/mol] molecular weight of xylene
        const Scalar M_a = 1e3*Air::molarMass(); // [g/mol] molecular weight of air
        const Scalar Tb_x = 412.0;        // [K] boiling temperature of xylene
        const Scalar sigma_a = 3.711;     // charact. length of air
        const Scalar T_scal_a = 78.6;     // [K] (molec. energy of attraction/Boltzmann constant)
        const Scalar V_B_x = 140.4;       // [cm^3/mol] LeBas molal volume of xylene
        const Scalar sigma_x = 1.18*pow(V_B_x, 0.333);     // charact. length of xylene
        const Scalar sigma_ax = 0.5*(sigma_a + sigma_x);
        const Scalar T_scal_x = 1.15*Tb_x;
        const Scalar T_scal_ax = sqrt(T_scal_a*T_scal_x);

        using std::max;
        Scalar T_star = temperature/T_scal_ax;
        T_star = max(T_star, 1e-5); // regularization

        const Scalar Omega = 1.06036/pow(T_star, 0.1561) + 0.193/exp(T_star*0.47635)
            + 1.03587/exp(T_star*1.52996) + 1.76474/exp(T_star*3.89411);
        const Scalar B_ = 0.00217 - 0.0005*sqrt(1.0/M_a + 1.0/M_x);
        const Scalar Mr = (M_a + M_x)/(M_a*M_x);
        const Scalar D_ax = (B_*pow(temperature,1.5)*sqrt(Mr))
                           /(1e-5*pressure*power(sigma_ax, 2)*Omega); // [cm^2/s]

        return D_ax*1e-4;   //  [m^2/s]
    }

    /*!
     * \brief Diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for air and xylene in liquid water.
     * \param temperature temperature in \f$\mathrm{[K]}\f$
     * \param pressure pressure in \f$\mathrm{[Pa]}\f$
     *
     * \note Returns just an order of magnitude.
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        return 1e-9;
    }
};

} // end namespace Dumux::BinaryCoeff

#endif
