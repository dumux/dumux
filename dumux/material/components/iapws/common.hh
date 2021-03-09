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
 * \ingroup IAPWS
 * \brief Implements relations common for all regions of the IAPWS '97
 *        formulation.
 * See:
 *
 * IAPWS: "Revised Release on the IAPWS Industrial Formulation
 * 1997 for the Thermodynamic Properties of Water and Steam",
 * http://www.iapws.org/relguide/IF97-Rev.pdf
 */
#ifndef DUMUX_IAPWS_COMMON_HH
#define DUMUX_IAPWS_COMMON_HH

#include <cmath>
#include <iostream>

#include <dune/common/math.hh>

#include <dumux/material/constants.hh>

namespace Dumux::IAPWS {

/*!
 * \ingroup IAPWS
 * \brief Implements relations which are common for all regions of the IAPWS '97
 *        formulation.
 *
 * \tparam Scalar The type used for scalar values
 *
 * See:
 *
 * IAPWS: "Revised Release on the IAPWS Industrial Formulation
 * 1997 for the Thermodynamic Properties of Water and Steam",
 * http://www.iapws.org/relguide/IF97-Rev.pdf
 */
template <class Scalar>
class Common
{
public:
    //! The molar mass of water \f$\mathrm{[kg/mol]}\f$
    static constexpr Scalar molarMass = 18.01518e-3;

    //! Specific gas constant of water \f$\mathrm{[J/(kg*K)]}\f$
    static constexpr Scalar Rs = Constants<Scalar>::R / molarMass;

    //! Critical temperature of water \f$\mathrm{[K]}\f$
    static constexpr Scalar criticalTemperature = 647.096;

    //! Critical pressure of water \f$\mathrm{[Pa]}\f$
    static constexpr Scalar criticalPressure = 22.064e6;

    //! Critical molar volume of water \f$\mathrm{[m^3/mol]}\f$
    static constexpr Scalar criticalMolarVolume = molarMass / 322.0;

    //! The acentric factor of water \f$\mathrm{[-]}\f$
    static constexpr Scalar acentricFactor = 0.344;

    //! Density of water at the critical point \f$\mathrm{[kg/m^3]}\f$
    static constexpr Scalar criticalDensity = 322;

    //! Triple temperature of water \f$\mathrm{[K]}\f$
    static constexpr Scalar tripleTemperature = 273.16;

    //! Triple pressure of water \f$\mathrm{[Pa]}\f$
    static constexpr Scalar triplePressure = 611.657;

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure water.
     *
     * This relation is valid for all regions of the IAPWS '97
     * formulation.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param rho density of component in \f$\mathrm{[kg/m^3]}\f$
     *
     * See:
     *
     * IAPWS: "Release on the IAPWS Formulation 2008 for the Viscosity
     * of Ordinary Water Substance", http://www.iapws.org/relguide/visc.pdf
     */
    static Scalar viscosity(Scalar temperature, Scalar rho)
    {
        const Scalar rhoBar = rho/322.0;
        const Scalar TBar = temperature/criticalTemperature;

        // muBar = muBar_1
        constexpr Scalar Hij[6][7] = {
            { 5.20094e-1, 2.22531e-1,-2.81378e-1, 1.61913e-1,-3.25372e-2, 0, 0 },
            { 8.50895e-2, 9.99115e-1,-9.06851e-1, 2.57399e-1, 0, 0, 0 },
            {-1.08374 , 1.88797 ,-7.72479e-1, 0, 0, 0, 0 },
            {-2.89555e-1, 1.26613 ,-4.89837e-1, 0, 6.98452e-2, 0,-4.35673e-3 },
            {          0, 0,-2.57040e-1, 0, 0, 8.72102e-3, 0 },
            {          0, 1.20573e-1, 0, 0, 0, 0,-5.93264e-4 }
        };

        Scalar tmp, tmp2, tmp3 = 1;
        Scalar muBar = 0;
        for (int i = 0; i <= 5; ++i) {
            tmp = 0;
            tmp2 = 1;
            for (int j = 0; j <= 6; ++j) {
                tmp += Hij[i][j]*tmp2;
                tmp2 *= (rhoBar - 1);
            }
            muBar += tmp3 * tmp;
            tmp3 *= 1.0/TBar - 1;
        }
        using std::exp;
        muBar *= rhoBar;
        muBar = exp(muBar);

        // muBar *= muBar_0
        using std::sqrt;
        muBar  *= 100*sqrt(TBar);
        constexpr Scalar H[4] = {
            1.67752, 2.20462, 0.6366564, -0.241605
        };

        tmp = 0, tmp2 = 1;
        for (int i = 0; i < 4; ++i) {
            tmp += H[i]/tmp2;
            tmp2 *= TBar;
        }
        muBar /= tmp;

        return 1e-6*muBar;
    }

    /*!
     * \brief Thermal conductivity \f$\mathrm{[[W/(m*K)]}\f$ water (IAPWS) .
     *
     * Implementation taken from:
     * freesteam - IAPWS-IF97 steam tables library
     * copyright (C) 2004-2009  John Pye
     *
     * Appendix B: Recommended Interpolating equation for Industrial Use
     * see http://www.iapws.org/relguide/thcond.pdf
     *
     * \param T absolute temperature in \f$\mathrm{[K]}\f$
     * \param rho density of water in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar thermalConductivityIAPWS(const Scalar T, const Scalar rho)
    {
        constexpr Scalar thcond_tstar   = 647.26 ;
        constexpr Scalar thcond_rhostar = 317.7 ;
        /*static constexpr Scalar thcond_kstar   = 1.0 ;*/

        constexpr Scalar thcond_b0      = -0.397070 ;
        constexpr Scalar thcond_b1      = 0.400302 ;
        constexpr Scalar thcond_b2      = 1.060000 ;
        constexpr Scalar thcond_B1      = -0.171587 ;
        constexpr Scalar thcond_B2      = 2.392190 ;

        constexpr Scalar thcond_c1      = 0.642857 ;
        constexpr Scalar thcond_c2      = -4.11717 ;
        constexpr Scalar thcond_c3      = -6.17937 ;
        constexpr Scalar thcond_c4      = 0.00308976 ;
        constexpr Scalar thcond_c5      = 0.0822994 ;
        constexpr Scalar thcond_c6      = 10.0932 ;

        constexpr Scalar thcond_d1      = 0.0701309 ;
        constexpr Scalar thcond_d2      = 0.0118520 ;
        constexpr Scalar thcond_d3      = 0.00169937 ;
        constexpr Scalar thcond_d4      = -1.0200 ;
        constexpr unsigned int thcond_a_count = 4;
        constexpr Scalar thcond_a[thcond_a_count] = {
            0.0102811
            ,0.0299621
            ,0.0156146
            ,-0.00422464
        };

        const Scalar Tbar = T / thcond_tstar;
        const Scalar rhobar = rho / thcond_rhostar;

        /* fast implementation... minimised calls to 'pow' routine... */
        using std::sqrt;
        const Scalar Troot = sqrt(Tbar);
        Scalar Tpow = Troot;
        Scalar lam = 0;

        for(unsigned int k = 0; k < thcond_a_count; ++k) {
            lam += thcond_a[k] * Tpow;
            Tpow *= Tbar;
        }

        using std::exp;
        lam += thcond_b0 + thcond_b1
                * rhobar + thcond_b2
                * exp(thcond_B1 * ((rhobar + thcond_B2)*(rhobar + thcond_B2)));

        using std::abs;
        using std::pow;
        using Dune::power;
        const Scalar DTbar = abs(Tbar - 1) + thcond_c4;
        const Scalar DTbarpow = pow(DTbar, 3./5);
        const Scalar Q = 2. + thcond_c5 / DTbarpow;
        const Scalar S = (Tbar >= 1) ? 1. / DTbar : thcond_c6 / DTbarpow;

        const Scalar rhobar18 = pow(rhobar, 1.8);
        const Scalar rhobarQ = pow(rhobar, Q);

        lam +=
            (thcond_d1 / power(Tbar,10) + thcond_d2) * rhobar18 *
                exp(thcond_c1 * (1 - rhobar * rhobar18))
            + thcond_d3 * S * rhobarQ *
                exp((Q/(1+Q))*(1 - rhobar*rhobarQ))
            + thcond_d4 *
                exp(thcond_c2 * power(Troot,3) + thcond_c3 / power(rhobar,5));
        return /*thcond_kstar * */ lam;
    }
};

} // end namespace Dumux::IAPWS

#endif
