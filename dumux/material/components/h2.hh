// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
#include <dumux/material/components/shomate.hh>

namespace Dumux::Components {

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
    using ShomateMethod = Dumux::ShomateMethod<Scalar, 3>; // three regions

public:
    static const ShomateMethod shomateMethod;

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
     *        Shomate Equation is used for a temperature range of 298K to 6000K.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasEnthalpy(Scalar temperature,
                                    Scalar pressure)
    {
        const auto h = shomateMethod.enthalpy(temperature); // KJ/mol
        return h * 1e3 / molarMass(); // J/kg
    }

    /*!
     * \brief Specific isobaric heat capacity \f$\mathrm{[J/(kg*K)]}\f$ of pure
     *        hydrogen gas.
     *        Shomate Equation is used for a temperature range of 298K to 6000K.
     *
     * \param T temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See: R. Reid, et al. (1987, pp 154, 657, 665) \cite reid1987
     */
    static const Scalar gasHeatCapacity(Scalar T,
                                        Scalar pressure)
    {
        const auto cp = shomateMethod.heatCapacity(T); // J/(mol K)
        return cp / molarMass(); // J/(kg K)
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

        // conversion from micro poise to Pa s
        return mu/1e6 / 10;
    }
};

/*!
 * \brief Shomate parameters for hydrogen published by NIST  \cite NIST
 * https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Units=SI&Mask=1&Type=JANAFG&Table=on#JANAFG
 * First row defines the temperature ranges, further rows give the parameters (A,B,C,D,E,F,G,H) for the respective temperature ranges.
 */
template <class Scalar>
const typename H2<Scalar>::ShomateMethod H2<Scalar>::shomateMethod{
    /*temperature*/{298.0, 1000.0, 2500.0, 6000.0},
    typename H2<Scalar>::ShomateMethod::Coefficients{{
        {33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974, 0.0},
        {18.563083, 12.257357, -2.859786, 0.268238, 1.97799, -1.147438, 156.288133, 0.0},
        {43.41356, -4.293079, 1.272428, -0.096876, -20.533862, -38.515158, 162.081354, 0.0}
    }}
};

} // end namespace Dumux::Components

#endif
