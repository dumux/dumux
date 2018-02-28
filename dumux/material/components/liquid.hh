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
 * \brief Interface for components that have a liquid state
 */
#ifndef DUMUX_COMPONENT_LIQUID_HH
#define DUMUX_COMPONENT_LIQUID_HH

namespace Dumux {
namespace Components {

template<class Scalar, class Component>
class Liquid
{
public:
    /*!
     * \brief the component has a liquid state if it derives from Liquid
     */
    static constexpr bool hasLiquidState()
    { return true; }

    /*!
     * \brief Returns true if the liquid phase is assumed to be compressible
     */
    static constexpr bool liquidIsCompressible()
    { return Component::liquidIsCompressible(); }

    /*!
     * \brief Returns true if the liquid phase viscostiy is constant
     */
    static constexpr bool liquidViscosityIsConstant()
    { return Component::liquidViscosityIsConstant(); }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the liquid component at a given pressure in
     *          \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    { return Component::liquidDensity(temperature, pressure); }

    /*!
     * \brief The dynamic liquid viscosity \f$\mathrm{[Pa*s]}\f$ of the pure component.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    { return Component::liquidViscosity(temperature, pressure); }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of the pure component in liquid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::liquidEnthalpy()"); }

    /*!
     * \brief Specific internal energy \f$\mathrm{[J/kg]}\f$ of pure the pure component in liquid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar liquidInternalEnergy(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::liquidInternalEnergy()"); }

    /*!
     * \brief Thermal conductivity of the component \f$\mathrm{[W/(m*K)]}\f$ as a liquid.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidThermalConductivity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::liquidThermalConductivity()"); }

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a liquid.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidHeatCapacity(Scalar temperature, Scalar pressure)
    { DUNE_THROW(Dune::NotImplemented, "Component::liquidHeatCapacity()"); }
};

} // end namespace Components
} // end namespace Dumux

#endif
