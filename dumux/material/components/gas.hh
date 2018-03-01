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
 * \ingroup Components
 * \brief Interface for components that have a gas state
 */
#ifndef DUMUX_COMPONENT_GAS_HH
#define DUMUX_COMPONENT_GAS_HH

namespace Dumux {
namespace Components {

template<class Scalar, class Component>
class Gas
{
public:
    /*!
     * \brief the component has a gas state if it derives from Gas
     */
    static constexpr bool hasGasState()
    { return true; }

    /*!
     * \brief Returns true if the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
    { return Component::gasIsCompressible(); }

    /*!
     * \brief Returns true if the gas phase viscostiy is constant
     */
    static constexpr bool gasViscosityIsConstant()
    { return Component::gasViscosityIsConstant(); }

    /*!
     * \brief Returns true if the gas phase is assumed to be ideal
     */
    static constexpr bool gasIsIdeal()
    { return Component::gasIsCompressible(); }

    /*!
     * \brief The density in \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure in
     *          \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::gasDensity()"); }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of the pure component in gas.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::gasEnthalpy()"); }

    /*!
     * \brief Specific internal energy \f$\mathrm{[J/kg]}\f$ of the pure component in gas.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasInternalEnergy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::gasInternalEnergy()"); }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of the pure component at a given pressure in
     *          \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::gasViscosity()"); }

    /*!
     * \brief Thermal conductivity of the component \f$\mathrm{[W/(m*K)]}\f$ as a gas.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasThermalConductivity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::gasThermalConductivity()"); }

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a gas.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasHeatCapacity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::gasHeatCapacity()"); }

};

} // end namespace Components
} // end namespace Dumux

#endif
