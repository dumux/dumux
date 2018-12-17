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
 * \ingroup Components
 * \brief Properties of pure molecular hydrogen \f$H_2\f$.
 */
#ifndef DUMUX_H2_HH
#define DUMUX_H2_HH

#include <dumux/material/idealgas.hh>

#include <cmath>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/gas.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief Properties of pure molecular hydrogen \f$H_2\f$.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class H2
: public Components::Base<Scalar, H2<Scalar> >
, public Components::Gas<Scalar, H2<Scalar> >
{
    using IdealGas = Dumux::IdealGas<Scalar>;

public:
    /*!
     * \brief A human readable name for the \f$H_2\f$.
     */
    static std::string name()
    { return "H2"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of molecular hydrogen.
     */
    static constexpr Scalar molarMass()
    { return 2.01588e-3; }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of molecular hydrogen.
     */
    static Scalar criticalTemperature()
    { return 33.2; /* [K] */ }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of molecular hydrogen.
     */
    static Scalar criticalPressure()
    { return 13.0e5; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at molecular hydrogen's triple point.
     */
    static Scalar tripleTemperature()
    { return 14.0; /* [K] */ }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure molecular hydrogen
     *        at a given temperature.
     *
     *\param temperature temperature of component in \f$\mathrm{[K]}\f$
     *
     * Taken from:
     *
     * See: R. Reid, et al. (1987, pp 208-209, 669) \cite reid1987
     *
     * \todo implement the Gomez-Thodos approach...
     */
    static Scalar vaporPressure(Scalar temperature)
    {
        if (temperature > criticalTemperature())
            return criticalPressure();
        if (temperature < tripleTemperature())
            return 0; // H2 is solid: We don't take sublimation into
                      // account

        // antoine equatuion
        const Scalar A = -7.76451;
        const Scalar B = 1.45838;
        const Scalar C = -2.77580;

        using std::exp;
        return 1e5 * exp(A - B/(temperature + C));
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of \f$H_2\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        // Assume an ideal gas
        return IdealGas::density(molarMass(), temperature, pressure);
    }

    /*!
     * \brief The molar density of \f$H_2\f$ in \f$\mathrm{[mol/m^3]}\f$,
     *   depending on pressure and temperature.
     * \param temperature The temperature of the gas
     * \param pressure The pressure of the gas
     */
    static Scalar gasMolarDensity(Scalar temperature, Scalar pressure)
    { return IdealGas::molarDensity(temperature, pressure); }

    /*!
     * \brief Returns true if the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true if the gas phase is assumed to be ideal
     */
    static constexpr bool gasIsIdeal()
    { return true; }

    /*!
     * \brief The pressure of gaseous \f$H_2\f$ in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        // Assume an ideal gas
        return IdealGas::pressure(temperature, density/molarMass());
    }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of pure hydrogen gas.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasEnthalpy(Scalar temperature,
                                    Scalar pressure)
    {
        return gasHeatCapacity(temperature, pressure) * temperature;
    }

    /*!
     * \brief Specific isobaric heat capacity \f$\mathrm{[J/(kg*K)]}\f$ of pure
     *        hydrogen gas.
     *
     * This is equivalent to the partial derivative of the specific
     * enthalpy to the temperature.
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See: R. Reid, et al. (1987, pp 154, 657, 665) \cite reid1987
     */
    static const Scalar gasHeatCapacity(Scalar T,
                                        Scalar pressure)
    {
        // method of Joback
        const Scalar cpVapA = 27.14;
        const Scalar cpVapB = 9.273e-3;
        const Scalar cpVapC = -1.381e-5;
        const Scalar cpVapD = 7.645e-9;

        return
            1/molarMass()* // conversion from [J/(mol*K)] to [J/(kg*K)]
            (cpVapA + T*
              (cpVapB/2 + T*
                (cpVapC/3 + T*
                 (cpVapD/4))));
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of \f$H_2\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * See: R. Reid, et al.: The Properties of Gases and Liquids,
     * 4th edition (1987, pp 396-397, 667) \cite reid1987 <BR>
     * 5th edition (2001, pp 9.7-9.8 (omega and V_c taken from p. A.19)) \cite poling2001
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        const Scalar Tc = criticalTemperature();
        const Scalar Vc = 65.0; // critical specific volume [cm^3/mol]
        const Scalar omega = -0.216; // accentric factor
        const Scalar M = molarMass() * 1e3; // molar mas [g/mol]
        const Scalar dipole = 0.0; // dipole moment [debye]

        using std::sqrt;
        Scalar mu_r4 = 131.3 * dipole / sqrt(Vc * Tc);
        mu_r4 *= mu_r4;
        mu_r4 *= mu_r4;

        using std::pow;
        using std::exp;
        Scalar Fc = 1 - 0.2756*omega + 0.059035*mu_r4;
        Scalar Tstar = 1.2593 * temperature/Tc;
        Scalar Omega_v =
            1.16145*pow(Tstar, -0.14874) +
            0.52487*exp(- 0.77320*Tstar) +
            2.16178*exp(- 2.43787*Tstar);
        Scalar mu = 40.785*Fc*sqrt(M*temperature)/(pow(Vc, 2./3)*Omega_v);

        // convertion from micro poise to Pa s
        return mu/1e6 / 10;
    }
};

} // end namespace Components

} // end namespace Dumux

#endif
