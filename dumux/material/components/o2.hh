// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief Properties of pure molecular oxygen \f$O_2\f$.
 */
#ifndef DUMUX_O2_HH
#define DUMUX_O2_HH

#include <dumux/material/idealgas.hh>

#include <cmath>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/gas.hh>
#include <dumux/material/components/shomate.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief Properties of pure molecular oxygen \f$O_2\f$.
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class O2
: public Components::Base<Scalar, O2<Scalar> >
, public Components::Gas<Scalar, O2<Scalar> >
{
    using IdealGas = Dumux::IdealGas<Scalar>;
    using ShomateMethod = Dumux::ShomateMethod<Scalar, 3>; // three regions

public:
    static const ShomateMethod shomateMethod;

    /*!
     * \brief A human readable name for the \f$O_2\f$.
     */
    static std::string name()
    { return "O2"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of molecular oxygen.
     */
    static constexpr Scalar molarMass()
    { return 32e-3; }

    /*!
     * \brief Returns the critical temperature in \f$\mathrm{[K]}\f$ of molecular oxygen.
     */
    static constexpr Scalar criticalTemperature()
    { return 154.581; /* [K] */ }

    /*!
     * \brief Returns the critical pressure in \f$\mathrm{[Pa]}\f$ of molecular oxygen.
     */
    static constexpr Scalar criticalPressure()
    { return 5.0804e6; /* [N/m^2] */ }

    /*!
     * \brief Returns the temperature in \f$\mathrm{[K]}\f$ at molecular oxygen's triple point.
     */
    static constexpr Scalar tripleTemperature()
    { return 54.359; /* [K] */ }

    /*!
     * \brief Returns the pressure in \f$\mathrm{[Pa]}\f$ at molecular oxygen's triple point.
     */
    static constexpr Scalar triplePressure()
    { return 148.0; /* [N/m^2] */ }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure molecular oxygen
     *        at a given temperature.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     *
     * Taken from:
     *
     * R. Prydz (1972, pp. 1-4) \cite prydz1972
     */
    static Scalar vaporPressure(Scalar T)
    {
        if (T > criticalTemperature())
            return criticalPressure();
        if (T < tripleTemperature())
            return 0; // O2 is solid: We don't take sublimation into account

        // vapor pressure between tripe and critical points.  See the
        // paper of Prydz for a discussion
        Scalar X =
            (1 - tripleTemperature()/T) /
            (1 - tripleTemperature()/criticalTemperature());
        const Scalar A = 7.568956;
        const Scalar B = 5.004836;
        const Scalar C = -2.137460;
        const Scalar D = 3.454481;
        const Scalar epsilon = 1.514;

        using std::exp;
        using std::pow;
        return triplePressure()*exp(X*(A + X*(B + C*X) + D*pow(1 - X, epsilon)));
    }

    /*!
     * \brief Returns true if the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
    { return true; }

    /*!
     * \brief The density in \f$\mathrm{[kg/m^3]}\f$ of pure \f$O_2\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * \todo: density liquid oxygen
     */
    static constexpr Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        // Assume an ideal gas
        return IdealGas::density(molarMass(), temperature, pressure);
    }

    /*!
     * \brief The molar density of pure \f$O_2\f$ in \f$\mathrm{[mol/m^3]}\f$,
     *   depending on pressure and temperature.
     * \param temperature The temperature of the gas
     * \param pressure The pressure of the gas
     */
    static Scalar gasMolarDensity(Scalar temperature, Scalar pressure)
    { return IdealGas::molarDensity(temperature, pressure); }

    /*!
     * \brief Returns true if the gas phase is assumed to be ideal
     */
    static constexpr bool gasIsIdeal()
    { return true; }

    /*!
     * \brief The pressure of gaseous \f$O_2\f$ in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static constexpr Scalar gasPressure(Scalar temperature, Scalar density)
    {
        // Assume an ideal gas
        return IdealGas::pressure(temperature, density/molarMass());
    }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of pure oxygen gas.
     * Shomate Equation is used for a temperature range of 100K to 6000K.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasEnthalpy(Scalar temperature,
                              Scalar pressure)
    {
        const auto h = shomateMethod.enthalpy(temperature); // KJ/mol
        return h * 1e3 / molarMass(); // J/kg
    }

    /*!
     * \brief Specific isobaric heat capacity \f$\mathrm{[J/(kg*K)]}\f$ of pure oxygen gas.
     * Shomate Equation is used for a temperature range of 100K to 6000K.
     *
     * \param T absolute temperature in \f$\mathrm{[K]}\f$
     * \param pressure of the phase in \f$\mathrm{[Pa]}\f$
     *
     * See: R. Reid, et al. (1987, pp 154, 657, 665) \cite reid1987
     */
    static Scalar gasHeatCapacity(Scalar T,
                                  Scalar pressure)
    {
        const auto cp = shomateMethod.heatCapacity(T); // J/(mol K)
        return cp / molarMass(); // J/(kg K)
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of \f$O_2\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * See: R. Reid, et al. (1987, pp 396-397, 664) \cite reid1987
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        const Scalar Tc = criticalTemperature();
        const Scalar Vc = 73.4; // critical specific volume [cm^3/mol]
        const Scalar omega = 0.025; // accentric factor
        const Scalar M = molarMass() * 1e3; // molar mas [g/mol]
        const Scalar dipole = 0.0; // dipole moment [debye]

        using std::sqrt;
        Scalar mu_r4 = 131.3 * dipole / sqrt(Vc * Tc);
        mu_r4 *= mu_r4;
        mu_r4 *= mu_r4;

        Scalar Fc = 1 - 0.2756*omega + 0.059035*mu_r4;
        Scalar Tstar = 1.2593 * temperature/Tc;

        using std::pow;
        using std::exp;
        Scalar Omega_v =
            1.16145*pow(Tstar, -0.14874) +
            0.52487*exp(- 0.77320*Tstar) +
            2.16178*exp(- 2.43787*Tstar);
        Scalar mu = 40.785*Fc*sqrt(M*temperature)/(pow(Vc, 2./3)*Omega_v);

        // conversion from micro poise to Pa s
        return mu/1e6 / 10;
    }

    /*!
     * \brief Thermal conductivity \f$\mathrm{[[W/(m*K)]}\f$ of nitrogen.
     *
     * Isobaric Properties for Nitrogen and Oxygen in: NIST Standard
     * Reference Database Number 69, Eds. P.J. Linstrom and
     * W.G. Mallard evaluated at p=.1 MPa, does not
     * change dramatically with p and can be interpolated linearly with temperature
     *
     * \param temperature absolute temperature in \f$\mathrm{[K]}\f$
     * \param pressure of the phase in \f$\mathrm{[Pa]}\f$
     */
    static constexpr Scalar gasThermalConductivity(Scalar temperature, Scalar pressure)
    {
        return 8.044e-5 * (temperature - 273.15) + 0.024486;
    }
};

/*!
 * \brief Shomate parameters for oxygen published by NIST  \cite NIST
 * https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Units=SI&Mask=1&Type=JANAFG&Table=on#JANAFG
 * First row defines the temperature ranges, further rows give the parameters (A,B,C,D,E,F,G,H) for the respective temperature ranges.
 */
template <class Scalar>
const typename O2<Scalar>::ShomateMethod O2<Scalar>::shomateMethod{
    /*temperature*/{100.0, 700.0, 2000.0, 6000.0},
    typename O2<Scalar>::ShomateMethod::Coefficients{{
        {31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 246.7945, 0.0},
        {30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 236.1663, 0.0},
        {20.91111, 10.72071, -2.020498, 0.146449, 9.245722, 5.337651, 237.6185, 0.0}
    }}
};

} // end namespace Components
} // end namespace Dumux

#endif
