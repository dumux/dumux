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
 * \brief Material properties of pure water \f$H_2O\f$.
 */
#ifndef DUMUX_H2O_HH
#define DUMUX_H2O_HH

#include <cmath>
#include <cassert>

#include <dumux/material/idealgas.hh>
#include <dumux/common/exceptions.hh>

#include "iapws/common.hh"
#include "iapws/region1.hh"
#include "iapws/region2.hh"
#include "iapws/region4.hh"

#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/gas.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief Material properties of pure water \f$H_2O\f$.
 * \tparam Scalar The type used for scalar values
 *
 * See:
 * IAPWS: "Revised Release on the IAPWS Industrial Formulation
 * 1997 for the Thermodynamic Properties of Water and Steam",
 * http://www.iapws.org/relguide/IF97-Rev.pdf \cite IAPWS1997
 */
template <class Scalar>
class H2O
: public Components::Base<Scalar, H2O<Scalar> >
, public Components::Liquid<Scalar, H2O<Scalar> >
, public Components::Gas<Scalar, H2O<Scalar> >
{

    using Common = IAPWS::Common<Scalar>;
    using Region1 = IAPWS::Region1<Scalar>;
    using Region2 = IAPWS::Region2<Scalar>;
    using Region4 = IAPWS::Region4<Scalar>;

    // specific gas constant of water
    static constexpr Scalar Rs = Common::Rs;

public:
    /*!
     * \brief A human readable name for the water.
     */
    static std::string name()
    { return "H2O"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of water.
     */
    static constexpr Scalar molarMass()
    { return Common::molarMass; }

    /*!
     * \brief The acentric factor \f$\mathrm{[-]}\f$ of water.
     */
    static constexpr Scalar acentricFactor()
    { return Common::acentricFactor; }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of water
     */
    static constexpr Scalar criticalTemperature()
    { return Common::criticalTemperature; }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of water.
     */
    static constexpr Scalar criticalPressure()
    { return Common::criticalPressure; }

    /*!
     * \brief Returns the molar volume \f$\mathrm{[m^3/mol]}\f$ of water at the critical point
     */
    static constexpr Scalar criticalMolarVolume()
    { return Common::criticalMolarVolume; }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at water's triple point.
     */
    static constexpr Scalar tripleTemperature()
    { return Common::tripleTemperature; }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at water's triple point.
     */
    static constexpr Scalar triplePressure()
    { return Common::triplePressure; }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of pure water
     *        at a given temperature.
     *
     *\param T temperature of component in \f$\mathrm{[K]}\f$
     *
     * See:
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf \cite IAPWS1997
     */
    static Scalar vaporPressure(Scalar T)
    {
        using std::min;
        using std::max;
        T = min(T, criticalTemperature());
        T = max(T,tripleTemperature());

        return Region4::saturationPressure(T);
    }

    /*!
     * \brief The vapor temperature in \f$\mathrm{[K]}\f$ of pure water
     *        at a given pressure.
     *
     *\param pressure pressure in \f$\mathrm{[Pa]}\f$
     *
     * See:
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf \cite IAPWS1997
     */
    static Scalar vaporTemperature(Scalar pressure)
    {
        using std::min;
        using std::max;
        pressure = min(pressure, criticalPressure());
        pressure = max(pressure, triplePressure());

        return Region4::vaporTemperature(pressure);
    }

    /*!
     * \brief Specific enthalpy of water steam \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf \cite IAPWS1997
     */
    static const Scalar gasEnthalpy(Scalar temperature,
                                    Scalar pressure)
    {
        Region2::checkValidityRange(temperature, pressure, "Enthalpy");

        // regularization
        if (pressure < triplePressure() - 100) {
            // We assume an ideal gas for low pressures to avoid the
            // 0/0 for the gas enthalpy at very low pressures. The
            // enthalpy of an ideal gas does not exhibit any
            // dependence on pressure, so we can just return the
            // specific enthalpy at the point of regularization, i.e.
            // the triple pressure - 100Pa
            return enthalpyRegion2_(temperature, triplePressure() - 100);
        }
        Scalar pv = vaporPressure(temperature);
        if (pressure > pv) {
            // the pressure is too high, in this case we use the slope
            // of the enthalpy at the vapor pressure to regularize
            Scalar dh_dp =
                Rs*temperature*
                Region2::tau(temperature)*
                Region2::dPi_dp(pv)*
                Region2::ddGamma_dTaudPi(temperature, pv);

            return
                enthalpyRegion2_(temperature, pv) +
                (pressure - pv)*dh_dp;
        }

        return enthalpyRegion2_(temperature, pressure);
    }

    /*!
     * \brief Specific enthalpy of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf \cite IAPWS1997
     */
    static const Scalar liquidEnthalpy(Scalar temperature,
                                       Scalar pressure)
    {
        Region1::checkValidityRange(temperature, pressure, "Enthalpy");

        // regularization
        Scalar pv = vaporPressure(temperature);
        if (pressure < pv) {
            // the pressure is too low, in this case we use the slope
            // of the enthalpy at the vapor pressure to regularize
            Scalar dh_dp =
                Rs * temperature*
                Region1::tau(temperature)*
                Region1::dPi_dp(pv)*
                Region1::ddGamma_dTaudPi(temperature, pv);

            return
                enthalpyRegion1_(temperature, pv) +
                (pressure - pv)*dh_dp;
        }

        return enthalpyRegion1_(temperature, pressure);
    }

    /*!
     * \brief Specific isobaric heat capacity of water steam \f$\mathrm{[J/(kg*K)}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf \cite IAPWS1997
     */
    static const Scalar gasHeatCapacity(Scalar temperature,
                                        Scalar pressure)
    {
        Region2::checkValidityRange(temperature, pressure, "Heat capacity");

        // regularization
        if (pressure < triplePressure() - 100) {
            return heatCap_p_Region2_(temperature, triplePressure() - 100);
        }
        Scalar pv = vaporPressure(temperature);
        if (pressure > pv) {
            // the pressure is too high, in this case we use the heat
            // cap at the vapor pressure to regularize
            return heatCap_p_Region2_(temperature, pv);
        }
        return heatCap_p_Region2_(temperature, pressure);
    }

    /*!
     * \brief Specific isobaric heat capacity of liquid water \f$\mathrm{[J/(kg*K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf \cite IAPWS1997
     */
    static const Scalar liquidHeatCapacity(Scalar temperature,
                                           Scalar pressure)
    {
        Region1::checkValidityRange(temperature, pressure, "Heat capacity");

        // regularization
        Scalar pv = vaporPressure(temperature);
        if (pressure < pv) {
            // the pressure is too low, in this case we use the heat cap at the vapor pressure to regularize
            return heatCap_p_Region1_(temperature, pv);
        }

        return heatCap_p_Region1_(temperature, pressure);
    }

    /*!
     * \brief Specific internal energy of liquid water \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf \cite IAPWS1997
     */
    static const Scalar liquidInternalEnergy(Scalar temperature,
                                             Scalar pressure)
    {
        Region1::checkValidityRange(temperature, pressure, "Internal energy");

        // regularization
        Scalar pv = vaporPressure(temperature);
        if (pressure < pv) {
            // the pressure is too low, in this case we use the slope
            // of the internal energy at the vapor pressure to
            // regularize

            /*
            // calculate the partial derivative of the internal energy
            // to the pressure at the vapor pressure.
            Scalar tau = Region1::tau(temperature);
            Scalar dGamma_dPi = Region1::dGamma_dPi(temperature, pv);
            Scalar ddGamma_dTaudPi = Region1::ddGamma_dTaudPi(temperature, pv);
            Scalar ddGamma_ddPi = Region1::ddGamma_ddPi(temperature, pv);
            Scalar pi = Region1::pi(pv);
            Scalar dPi_dp = Region1::dPi_dp(pv);
            Scalar du_dp =
                Rs*temperature*
                (tau*dPi_dp*ddGamma_dTaudPi + dPi_dp*dPi_dp*dGamma_dPi + pi*dPi_dp*ddGamma_ddPi);
            */

            // use a straight line for extrapolation. use forward
            // differences to calculate the partial derivative to the
            // pressure at the vapor pressure
            static const Scalar eps = 1e-7;
            Scalar uv = internalEnergyRegion1_(temperature, pv);
            Scalar uvPEps = internalEnergyRegion1_(temperature, pv + eps);
            Scalar du_dp = (uvPEps - uv)/eps;
            return uv + du_dp*(pressure - pv);
        }

        return internalEnergyRegion1_(temperature, pressure);
    }

    /*!
     * \brief Specific internal energy of steam and water vapor \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf  \cite IAPWS1997
     */
    static Scalar gasInternalEnergy(Scalar temperature, Scalar pressure)
    {
        Region2::checkValidityRange(temperature, pressure, "Internal energy");

        // regularization
        if (pressure < triplePressure() - 100) {
            // We assume an ideal gas for low pressures to avoid the
            // 0/0 for the internal energy of gas at very low
            // pressures. The enthalpy of an ideal gas does not
            // exhibit any dependence on pressure, so we can just
            // return the specific enthalpy at the point of
            // regularization, i.e.  the triple pressure - 100Pa, and
            // subtract the work required to change the volume for an
            // ideal gas.
            return
                enthalpyRegion2_(temperature, triplePressure() - 100)
                -
                Rs*temperature; // = p*v   for an ideal gas!
        }
        Scalar pv = vaporPressure(temperature);
        if (pressure > pv) {
            // the pressure is too high, in this case we use the slope
            // of the internal energy at the vapor pressure to
            // regularize

            /*
            // calculate the partial derivative of the internal energy
            // to the pressure at the vapor pressure.
            Scalar tau = Region2::tau(temperature);
            Scalar dGamma_dPi = Region2::dGamma_dPi(temperature, pv);
            Scalar ddGamma_dTaudPi = Region2::ddGamma_dTaudPi(temperature, pv);
            Scalar ddGamma_ddPi = Region2::ddGamma_ddPi(temperature, pv);
            Scalar pi = Region2::pi(pv);
            Scalar dPi_dp = Region2::dPi_dp(pv);
            Scalar du_dp =
                Rs*temperature*
                (tau*dPi_dp*ddGamma_dTaudPi + dPi_dp*dPi_dp*dGamma_dPi + pi*dPi_dp*ddGamma_ddPi);

            // use a straight line for extrapolation
            Scalar uv = internalEnergyRegion2_(temperature, pv);
            return uv + du_dp*(pressure - pv);
            */

            // use a straight line for extrapolation. use backward
            // differences to calculate the partial derivative to the
            // pressure at the vapor pressure
            static const Scalar eps = 1e-7;
            Scalar uv = internalEnergyRegion2_(temperature, pv);
            Scalar uvMEps = internalEnergyRegion2_(temperature, pv - eps);
            Scalar du_dp = (uv - uvMEps)/eps;
            return uv + du_dp*(pressure - pv);
        }

        return internalEnergyRegion2_(temperature, pressure);
    }

    /*!
     * \brief Specific isochoric heat capacity of liquid water \f$\mathrm{[J/(m^3*K)]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf \cite IAPWS1997
     */
    static const Scalar liquidHeatCapacityConstVolume(Scalar temperature,
                                                      Scalar pressure)
    {
        Region1::checkValidityRange(temperature, pressure, "Heat capacity for a constant volume");

        // regularization
        Scalar pv = vaporPressure(temperature);
        if (pressure < pv) {
            // the pressure is too low, in this case we use the heat cap at the vapor pressure to regularize

            return heatCap_v_Region1_(temperature, pv);
        }

        return heatCap_v_Region1_(temperature, pressure);
    }

    /*!
     * \brief Specific isochoric heat capacity of steam and water vapor \f$\mathrm{[J/(kg*K)}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf \cite IAPWS1997
     */
    static Scalar gasHeatCapacityConstVolume(Scalar temperature, Scalar pressure)
    {
        Region2::checkValidityRange(temperature, pressure, "Heat capacity for a constant volume");

        // regularization
        if (pressure < triplePressure() - 100) {
            return
                heatCap_v_Region2_(temperature, triplePressure() - 100);
        }
        Scalar pv = vaporPressure(temperature);
        if (pressure > pv) {
            return heatCap_v_Region2_(temperature, pv);
        }

        return heatCap_v_Region2_(temperature, pressure);
    }

    /*!
     * \brief Returns true if the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true if the liquid phase is assumed to be compressible
     */
    static constexpr bool liquidIsCompressible()
    { return true; }

    /*!
     * \brief The density of steam in \f$\mathrm{[kg/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf \cite IAPWS1997
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        Region2::checkValidityRange(temperature, pressure, "Density");

        // regularization
        if (pressure < triplePressure() - 100) {
            // We assume an ideal gas for low pressures to avoid the
            // 0/0 for the internal energy and enthalpy.
            Scalar rho0IAPWS = 1.0/volumeRegion2_(temperature,
                                                  triplePressure() - 100);
            Scalar rho0Id = IdealGas<Scalar>::density(molarMass(),
                                                      temperature,
                                                      triplePressure() - 100);
            return
                rho0IAPWS/rho0Id *
                IdealGas<Scalar>::density(molarMass(),
                                          temperature,
                                          pressure);
        }
        Scalar pv = vaporPressure(temperature);
        if (pressure > pv) {
            // the pressure is too high, in this case we use the slope
            // of the density energy at the vapor pressure to
            // regularize

            // calculate the partial derivative of the specific volume
            // to the pressure at the vapor pressure.
            const Scalar eps = pv*1e-8;
            Scalar v0 = volumeRegion2_(temperature, pv);
            Scalar v1 = volumeRegion2_(temperature, pv + eps);
            Scalar dv_dp = (v1 - v0)/eps;
            /*
            Scalar pi = Region2::pi(pv);
            Scalar dp_dPi = Region2::dp_dPi(pv);
            Scalar dGamma_dPi = Region2::dGamma_dPi(temperature, pv);
            Scalar ddGamma_ddPi = Region2::ddGamma_ddPi(temperature, pv);

            Scalar RT = Rs*temperature;
            Scalar dv_dp =
                RT/(dp_dPi*pv)
                *
                (dGamma_dPi + pi*ddGamma_ddPi - v0*dp_dPi/RT);
            */

            // calculate the partial derivative of the density to the
            // pressure at vapor pressure
            Scalar drho_dp = - 1/(v0*v0)*dv_dp;

            // use a straight line for extrapolation
            return 1.0/v0 + (pressure - pv)*drho_dp;
        }

        return 1.0/volumeRegion2_(temperature, pressure);
    }

    /*!
     *  \brief The molar density of steam in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar gasMolarDensity(Scalar temperature, Scalar pressure)
    { return gasDensity(temperature, pressure)/molarMass(); }

    /*!
     * \brief Returns true if the gas phase is assumed to be ideal
     */
    static constexpr bool gasIsIdeal()
    { return false; }

    /*!
     * \brief The pressure of steam in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density of component in \f$\mathrm{[kg/m^3]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf \cite IAPWS1997
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        // We use the newton method for this. For the initial value we
        // assume steam to be an ideal gas
        Scalar pressure = IdealGas<Scalar>::pressure(temperature, density/molarMass());
        Scalar eps = pressure*1e-7;

        Scalar deltaP = pressure*2;
        using std::abs;
        for (int i = 0; i < 5 && abs(pressure*1e-9) < abs(deltaP); ++i)
        {
            const Scalar f = gasDensity(temperature, pressure) - density;

            Scalar df_dp;
            df_dp = gasDensity(temperature, pressure + eps);
            df_dp -= gasDensity(temperature, pressure - eps);
            df_dp /= 2*eps;

            deltaP = - f/df_dp;

            pressure += deltaP;
        }

        return pressure;
    }

    /*!
     * \brief The density of pure water in \f$\mathrm{[kg/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf \cite IAPWS1997
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        Region1::checkValidityRange(temperature, pressure, "Density");

        // regularization
        Scalar pv = vaporPressure(temperature);
        if (pressure < pv) {
            // the pressure is too low, in this case we use the slope
            // of the density at the vapor pressure to regularize

            // calculate the partial derivative of the specific volume
            // to the pressure at the vapor pressure.
            const Scalar eps = pv*1e-8;
            Scalar v0 = volumeRegion1_(temperature, pv);
            Scalar v1 = volumeRegion1_(temperature, pv + eps);
            Scalar dv_dp = (v1 - v0)/eps;

            /*
            Scalar v0 = volumeRegion1_(temperature, pv);
            Scalar pi = Region1::pi(pv);
            Scalar dp_dPi = Region1::dp_dPi(pv);
            Scalar dGamma_dPi = Region1::dGamma_dPi(temperature, pv);
            Scalar ddGamma_ddPi = Region1::ddGamma_ddPi(temperature, pv);

            Scalar RT = Rs*temperature;
            Scalar dv_dp =
                RT/(dp_dPi*pv)
                *
                (dGamma_dPi + pi*ddGamma_ddPi - v0*dp_dPi/RT);
            */

            // calculate the partial derivative of the density to the
            // pressure at vapor pressure
            Scalar drho_dp = - 1/(v0*v0)*dv_dp;

            // use a straight line for extrapolation
            return 1.0/v0 + (pressure - pv)*drho_dp;
        }

        return 1/volumeRegion1_(temperature, pressure);
    }

    /*!
     * \brief The molar density of water in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar liquidMolarDensity(Scalar temperature, Scalar pressure)
    { return liquidDensity(temperature, pressure)/molarMass(); }

    /*!
     * \brief The pressure of liquid water in \f$\mathrm{[Pa]}\f$ at a given density and
     *        temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     *
     * See:
     *
     * IAPWS: "Revised Release on the IAPWS Industrial Formulation
     * 1997 for the Thermodynamic Properties of Water and Steam",
     * http://www.iapws.org/relguide/IF97-Rev.pdf \cite IAPWS1997
     */
    static Scalar liquidPressure(Scalar temperature, Scalar density)
    {
        // We use the newton method for this. For the initial value we
        // assume the pressure to be 10% higher than the vapor
        // pressure
        Scalar pressure = 1.1*vaporPressure(temperature);
        Scalar eps = pressure*1e-7;

        Scalar deltaP = pressure*2;

        using std::abs;
        for (int i = 0; i < 5 && abs(pressure*1e-9) < abs(deltaP); ++i) {
            Scalar f = liquidDensity(temperature, pressure) - density;

            Scalar df_dp;
            df_dp = liquidDensity(temperature, pressure + eps);
            df_dp -= liquidDensity(temperature, pressure - eps);
            df_dp /= 2*eps;

            deltaP = - f/df_dp;

            pressure += deltaP;
        }

        return pressure;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of steam.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * This method is only valid if pressure is below or at the vapor
     * pressure of water.
     *
     * See:
     * IAPWS: "Release on the IAPWS Formulation 2008 for the Viscosity
     * of Ordinary Water Substance", http://www.iapws.org/relguide/visc.pdf \cite cooper2008
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        Region2::checkValidityRange(temperature, pressure, "Viscosity");

        Scalar rho = gasDensity(temperature, pressure);
        return Common::viscosity(temperature, rho);
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure water.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * See:
     * IAPWS: "Release on the IAPWS Formulation 2008 for the Viscosity
     * of Ordinary Water Substance", http://www.iapws.org/relguide/visc.pdf \cite cooper2008
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        Region1::checkValidityRange(temperature, pressure, "Viscosity");

        Scalar rho = liquidDensity(temperature, pressure);
        return Common::viscosity(temperature, rho);
    }

    /*!
     * \brief Thermal conductivity \f$\mathrm{[[W/(m*K)]}\f$ of water (IAPWS) .
     *
     * Implementation taken from:
     * freesteam - IAPWS-IF97 steam tables library
     * copyright (C) 2004-2009  John Pye
     *
     * See:
     * IAPWS: "Release on the IAPWS Formulation 2011 for the Thermal Conductivity
     * of Ordinary Water Substance", http://www.iapws.org/relguide/ThCond.pdf
     * \cite IAPWS_ThCond
     *
     * \param temperature absolute temperature in \f$\mathrm{[K]}\f$
     * \param pressure of the phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidThermalConductivity( Scalar temperature,  Scalar pressure)
    {
        // Thermal conductivity of water is empirically fit.
        // Evaluating that fitting-function outside the area of validity does not make sense.
        if ( !(   (pressure <= 400e6 && (273.15 <= temperature) && (temperature <= 398.15))
               || (pressure <= 200e6 && (398.15 <  temperature) && (temperature <=  523.15))
               || (pressure <= 150e6 && (523.15 <  temperature) && (temperature <=  673.15))
               || (pressure <= 100e6 && (673.15 <  temperature) && (temperature <=  1073.15)) ))
        {
            DUNE_THROW(NumericalProblem,
                       "Evaluating the IAPWS fit function for thermal conductivity outside range of applicability."
                       "(T=" << temperature << ", p=" << pressure << ")");
        }

        Scalar rho = liquidDensity(temperature, pressure);
        return Common::thermalConductivityIAPWS(temperature, rho);
    }

    /*!
     * \brief Thermal conductivity \f$\mathrm{[[W/(m*K)]}\f$ of steam (IAPWS) .
     *
     * Implementation taken from:
     * freesteam - IAPWS-IF97 steam tables library
     * copyright (C) 2004-2009  John Pye
     *
     * See:
     * IAPWS: "Release on the IAPWS Formulation 2011 for the Thermal Conductivity
     * of Ordinary Water Substance", http://www.iapws.org/relguide/ThCond.pdf
     * \cite IAPWS_ThCond
     *
     * \param temperature absolute temperature in \f$\mathrm{[K]}\f$
     * \param pressure of the phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasThermalConductivity(const Scalar temperature, const Scalar pressure)
    {
        // Thermal conductivity of water is empirically fit.
        // Evaluating that fitting-function outside the area of validity does not make sense.
        if ( !(   (pressure <= 400e6 && (273.15 <= temperature) && (temperature <= 398.15))
               || (pressure <= 200e6 && (398.15 <  temperature) && (temperature <= 523.15))
               || (pressure <= 150e6 && (523.15 <  temperature) && (temperature <= 673.15))
               || (pressure <= 100e6 && (673.15 <  temperature) && (temperature <= 1073.15)) ))
        {
            DUNE_THROW(NumericalProblem,
                       "Evaluating the IAPWS fit function for thermal conductivity outside range of applicability."
                       "(T=" << temperature << ", p=" << pressure << ")");
        }

        Scalar rho = gasDensity(temperature, pressure);
        return Common::thermalConductivityIAPWS(temperature, rho);
    }

private:
    // the unregularized specific enthalpy for liquid water
    static constexpr Scalar enthalpyRegion1_(Scalar temperature, Scalar pressure)
    {
        return
            Region1::tau(temperature) *
            Region1::dGamma_dTau(temperature, pressure) *
            Rs*temperature;
    }

    // the unregularized specific isobaric heat capacity
    static constexpr Scalar heatCap_p_Region1_(Scalar temperature, Scalar pressure)
    {
        return
            - Region1::tau(temperature) * Region1::tau(temperature) *
            Region1::ddGamma_ddTau(temperature, pressure) *
            Rs;
    }

    // the unregularized specific isochoric heat capacity
    static Scalar heatCap_v_Region1_(Scalar temperature, Scalar pressure)
    {
        Scalar tau = Region1::tau(temperature);
        Scalar num = Region1::dGamma_dPi(temperature, pressure) - tau * Region1::ddGamma_dTaudPi(temperature, pressure);
        Scalar diff = num * num / Region1::ddGamma_ddPi(temperature, pressure);

        return
            - tau * tau *
            Region1::ddGamma_ddTau(temperature, pressure) * Rs +
            diff;
    }

    // the unregularized specific internal energy for liquid water
    static constexpr Scalar internalEnergyRegion1_(Scalar temperature, Scalar pressure)
    {
        return
            Rs * temperature *
            ( Region1::tau(temperature)*Region1::dGamma_dTau(temperature, pressure) -
              Region1::pi(pressure)*Region1::dGamma_dPi(temperature, pressure));
    }

    // the unregularized specific volume for liquid water
    static constexpr Scalar volumeRegion1_(Scalar temperature, Scalar pressure)
    {
        return
            Region1::pi(pressure)*
            Region1::dGamma_dPi(temperature, pressure) *
            Rs * temperature / pressure;
    }

    // the unregularized specific enthalpy for steam
    static constexpr Scalar enthalpyRegion2_(Scalar temperature, Scalar pressure)
    {
        return
            Region2::tau(temperature) *
            Region2::dGamma_dTau(temperature, pressure) *
            Rs*temperature;
    }

    // the unregularized specific internal energy for steam
    static constexpr Scalar internalEnergyRegion2_(Scalar temperature, Scalar pressure)
    {
        return
            Rs * temperature *
            ( Region2::tau(temperature)*Region2::dGamma_dTau(temperature, pressure) -
              Region2::pi(pressure)*Region2::dGamma_dPi(temperature, pressure));
    }

    // the unregularized specific isobaric heat capacity
    static constexpr Scalar heatCap_p_Region2_(Scalar temperature, Scalar pressure)
    {
        return
            - Region2::tau(temperature) * Region2::tau(temperature) *
            Region2::ddGamma_ddTau(temperature, pressure) *
            Rs;
    }

    // the unregularized specific isochoric heat capacity
    static Scalar heatCap_v_Region2_(Scalar temperature, Scalar pressure)
    {
        Scalar tau = Region2::tau(temperature);
        Scalar pi = Region2::pi(pressure);
        Scalar num = 1 + pi * Region2::dGamma_dPi(temperature, pressure) + tau * pi * Region2::ddGamma_dTaudPi(temperature, pressure);
        Scalar diff = num * num / (1 - pi * pi * Region2::ddGamma_ddPi(temperature, pressure));
        return
            - tau * tau *
            Region2::ddGamma_ddTau(temperature, pressure) * Rs
            - diff;
    }

    // the unregularized specific volume for steam
    static constexpr Scalar volumeRegion2_(Scalar temperature, Scalar pressure)
    {
        return
            Region2::pi(pressure)*
            Region2::dGamma_dPi(temperature, pressure) *
            Rs * temperature / pressure;
    }
}; // end class

template <class Scalar>
struct IsAqueous<H2O<Scalar>> : public std::true_type {};

} // end namespace Components

} // end namespace Dumux

#endif
