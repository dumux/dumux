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
 * \brief Binary coefficients for water and xylene.
 */
#ifndef DUMUX_BINARY_COEFF_H2O_XYLENE_HH
#define DUMUX_BINARY_COEFF_H2O_XYLENE_HH

#include <algorithm>

#include <dune/common/math.hh>

#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/xylene.hh>

namespace Dumux::BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and xylene.
 */
class H2O_Xylene
{
public:
    /*!
     * \brief Henry coefficient \f$\mathrm{[Pa]}\f$  for xylene in liquid water.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     *
     * See:
     *
     *  Sander (1999) \cite sander1999
     */
    template <class Scalar>
    static Scalar henry(Scalar temperature)
    {
        // after Sander
        constexpr Scalar sanderH = 1.5e-1;    //[M/atm]
        //conversion to our Henry definition
        Scalar dumuxH = sanderH / 101.325; // has now [(mol/m^3)/Pa]
        dumuxH *= 18.02e-6;  //multiplied by molar volume of reference phase = water
        return 1.0/dumuxH; // [Pa]
    }

    /*!
     * \brief Binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular water and xylene.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the pressure \f$\mathrm{[Pa]}\f$
     *
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        using H2O = Dumux::Components::H2O<Scalar>;
        using Xylene = Dumux::Components::Xylene<Scalar>;

        using std::clamp;
        temperature = clamp(temperature, 1e-9, 500.0); // regularization
        pressure = clamp(pressure, 0.0, 1e8); // regularization

        using std::sqrt;
        using std::exp;
        using std::pow;
        using Dune::power;
        constexpr Scalar M_x = 1e3*Xylene::molarMass(); // [g/mol] molecular weight of xylene
        constexpr Scalar M_w = 1e3*H2O::molarMass(); // [g/mol] molecular weight of water
        constexpr Scalar Tb_x = 412.9;        // [K] boiling temperature of xylene
        constexpr Scalar Tb_w = 373.15;       // [K] boiling temperature of water (at p_atm)
        constexpr Scalar V_B_w = 18.0;                // [cm^3/mol] LeBas molal volume of water
        const Scalar  sigma_w = 1.18*pow(V_B_w, 0.333);     // charact. length of air
        constexpr Scalar  T_scal_w = 1.15*Tb_w;     // [K] (molec. energy of attraction/Boltzmann constant)
        constexpr Scalar V_B_x = 140.4;       // [cm^3/mol] LeBas molal volume of xylene
        const Scalar sigma_x = 1.18*pow(V_B_x, 0.333);     // charact. length of xylene
        const Scalar sigma_wx = 0.5*(sigma_w + sigma_x);
        constexpr Scalar T_scal_x = 1.15*Tb_x;
        const Scalar T_scal_wx = sqrt(T_scal_w*T_scal_x);

        using std::max;
        Scalar T_star = temperature/T_scal_wx;
        T_star = max(T_star, 1e-5); // regularization

        const Scalar Omega = 1.06036/pow(T_star, 0.1561) + 0.193/exp(T_star*0.47635)
            + 1.03587/exp(T_star*1.52996) + 1.76474/exp(T_star*3.89411);
        const Scalar  B_ = 0.00217 - 0.0005*sqrt(1.0/M_w + 1.0/M_x);
        const Scalar Mr = (M_w + M_x)/(M_w*M_x);
        const Scalar D_wx = (B_*pow(temperature,1.6)*sqrt(Mr))
                           /(1e-5*pressure*power(sigma_wx, 2)*Omega); // [cm^2/s]

        return D_wx*1e-4; // [m^2/s]
    }

    /*!
     * \brief Diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for xylene in liquid water.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the pressure \f$\mathrm{[Pa]}\f$
     *
     * \note Returns just an order of magnitude.
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        return 1.e-9;
    }
};

} // end namespace Dumux::BinaryCoeff

#endif
