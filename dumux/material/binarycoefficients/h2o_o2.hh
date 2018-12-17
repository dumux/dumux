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
 * \brief Binary coefficients for water and oxygen.
 */
#ifndef DUMUX_BINARY_COEFF_H2O_O2_HH
#define DUMUX_BINARY_COEFF_H2O_O2_HH

#include <dumux/material/binarycoefficients/henryiapws.hh>
#include <dumux/material/binarycoefficients/fullermethod.hh>

#include <dumux/material/components/o2.hh>
#include <dumux/material/components/h2o.hh>

namespace Dumux {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and oxygen.
 */
class H2O_O2
{
public:
  /*!
     * \brief Henry coefficient \f$\mathrm{[Pa]}\f$  for molecular oxygen in liquid water.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     */
    template <class Scalar>
    static Scalar henry(Scalar temperature)
    {
        const Scalar E = 2305.0674;
        const Scalar F = -11.3240;
        const Scalar G = 25.3224;
        const Scalar H = -15.6449;

        return henryIAPWS(E, F, G, H, temperature);
    }

    /*!
     * \brief Binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular water and oxygen.
     *
     * Uses fullerMethod to determine the diffusion of water in nitrogen.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        using H2O = Dumux::Components::H2O<Scalar>;
        using O2 = Components::O2<Scalar>;

        // atomic diffusion volumes
        const Scalar SigmaNu[2] = { 13.1 /* H2O */,  16.3 /* O2 */ };
        // molar masses [g/mol]
        const Scalar M[2] = { H2O::molarMass()*1e3, O2::molarMass()*1e3 };

        return fullerMethod(M, SigmaNu, temperature, pressure);
    }

    /*!
     * \brief Diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular oxygen in liquid water.
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
     *
     * R. Reid et al. (1987, pp. 599) \cite reid1987 <BR>
     *
     * R. Ferrell, D. Himmelblau (1967, pp. 111-115) \cite ferrell1967
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        const Scalar Texp = 273.15 + 25; // [K]
        const Scalar Dexp = 2.2e-9; // [m^2/s]
        return Dexp * temperature/Texp;
    }
};

} // end namespace BinaryCoeff
} // end namespace Dumux

#endif
