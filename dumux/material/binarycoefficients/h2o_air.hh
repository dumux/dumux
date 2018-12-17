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
 * \brief Binary coefficients for water and air.
 */
#ifndef DUMUX_BINARY_COEFF_H2O_AIR_HH
#define DUMUX_BINARY_COEFF_H2O_AIR_HH

#include <cmath>

namespace Dumux {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and air.
 */
class H2O_Air
{
public:
    /*!
     * \brief Henry coefficient \f$\mathrm{[Pa]}\f$  for air in liquid water.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     *
     * Henry coefficient See:
     * Stefan Finsterle (1993, page 33 Formula (2.9)) \cite finsterle1993 <BR>
     * (fitted to data from Tchobanoglous & Schroeder, 1985 \cite tchobanoglous1985 )
     */
    template <class Scalar>
    static Scalar henry(Scalar temperature)
    {
      using std::exp;
      Scalar r = (0.8942+1.47*exp(-0.04394*(temperature-273.15)))*1.E-10;

      return 1./r;
    }

    /*!
     * \brief Binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular water and air
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     * Vargaftik: Tables on the thermophysical properties of liquids and gases.
     * John Wiley & Sons, New York, 1975. \cite vargaftik1975 <BR>
     * Walker, Sabey, Hampton: Studies of heat transfer and water migration in soils.
     * Dep. of Agricultural and Chemical Engineering, Colorado State University,
     * Fort Collins, 1981. \cite walker1981
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        const Scalar Theta=1.8;
        const Scalar Daw=2.13e-5;  /* reference value */
        const Scalar pg0=1.e5;     /* reference pressure */
        const Scalar T0=273.15;    /* reference temperature */
        Scalar Dgaw;

        using std::pow;
        Dgaw=Daw*(pg0/pressure)*pow((temperature/T0),Theta);

        return Dgaw;
    }

    /*!
     * Lacking better data on water-air diffusion in liquids, we use at the
     * moment the diffusion coefficient of the air's main component nitrogen!!
     * \brief Diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular nitrogen in liquid water.
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     *
     * The empirical equations for estimating the diffusion coefficient in
     * infinite solution which are presented in Reid, 1987 all show a
     * linear dependency on temperature. We thus simply scale the
     * experimentally obtained diffusion coefficient of Ferrell and
     * Himmelblau by the temperature.
     *
     * See:
     * R. Reid et al. (1987, pp. 599) \cite reid1987 <BR>
     * R. Ferrell, D. Himmelblau (1967, pp. 111-115) \cite ferrell1967
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        const Scalar Texp = 273.15 + 25; // [K]
        const Scalar Dexp = 2.01e-9; // [m^2/s]
        return Dexp * temperature/Texp;
    }
};

} // end namespace BinaryCoeff
} // end namespace Dumux

#endif
