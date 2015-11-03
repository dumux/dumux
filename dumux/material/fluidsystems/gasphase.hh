// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 *
 * \brief A gaseous phase consisting of a single component.
 */
#ifndef DUMUX_GAS_PHASE_HH
#define DUMUX_GAS_PHASE_HH

namespace Dumux
{

/*!
 * \ingroup Fluidsystems
 * \brief gaseous phase consisting of a single component
 */
template <class Scalar, class ComponentT>
class GasPhase
{
public:
    typedef ComponentT Component;
    /*!
     * \brief A human readable name for the component.
     */
    static const char *name()
    { return Component::name(); }

    /*!
     * \brief Returs whether the fluid is a liquid
     */
    static bool isLiquid()
    { return false; }

    /*!
     * \brief Returns true iff the fluid is assumed to be compressible
     */
    static bool isCompressible()
    { return Component::gasIsCompressible(); }

    /*!
     * \brief Returns true iff the fluid is assumed to be an ideal gas
     */
    static bool isIdealGas()
    { return Component::gasIsIdeal(); }

    /*!
     * \brief The mass in \f$\mathrm{[kg]}\f$ of one mole of the component.
     */
    static Scalar molarMass()
    {  return Component::molarMass(); }

    /*!
     * \brief Returns the critical temperature in \f$\mathrm{[K]}\f$ of the component
     */
    static Scalar criticalTemperature()
    {  return Component::criticalTemperature(); }

    /*!
     * \brief Returns the critical pressure in \f$\mathrm{[Pa]}\f$ of the component
     */
    static Scalar criticalPressure()
    {  return Component::criticalPressure(); }

    /*!
     * \brief Returns the temperature in \f$\mathrm{[K]}\f$ at the component's triple point.
     */
    static Scalar tripleTemperature()
    {  return Component::tripleTemperature(); }

    /*!
     * \brief Returns the pressure in \f$\mathrm{[Pa]}\f$ at the component's triple point.
     */
    static Scalar triplePressure()
    { return Component::triplePressure(); }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of the component at a given
     *        temperature.
     * \param T temperature \f$\mathrm{[K]}\f$
     */
    static Scalar vaporPressure(Scalar T)
    { return Component::vaporPressure(T); }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure and temperature.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar density(Scalar temperature, Scalar pressure)
    {  return Component::gasDensity(temperature, pressure); }

    /*!
     * \brief The pressure \f$\mathrm{[Pa]}\f$ of the component at a given density and temperature.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param density The given density \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar pressure(Scalar temperature, Scalar density)
    {  return Component::gasPressure(temperature, density); }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of the pure component as a gas.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static const Scalar enthalpy(Scalar temperature, Scalar pressure)
    {  return Component::gasEnthalpy(temperature, pressure); }

    /*!
     * \brief Specific internal energy \f$\mathrm{[J/kg]}\f$ of the pure component as a gas.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static const Scalar internalEnergy(Scalar temperature, Scalar pressure)
    { return Component::gasInternalEnergy(temperature, pressure); }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the pure component at a given pressure and temperature.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar viscosity(Scalar temperature, Scalar pressure)
    {  return Component::gasViscosity(temperature, pressure); }

    /*!
     * \brief Thermal conductivity of the fluid \f$\mathrm{[W/(m K)]}\f$.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar thermalConductivity(Scalar temperature, Scalar pressure)
    { return Component::gasThermalConductivity(temperature, pressure); }

    /*!
     * \brief Specific isobaric heat capacity of the fluid \f$\mathrm{[J/(kg K)]}\f$.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar heatCapacity(Scalar temperature, Scalar pressure)
    { return Component::gasHeatCapacity(temperature, pressure); }
};
} // namespace

#endif
