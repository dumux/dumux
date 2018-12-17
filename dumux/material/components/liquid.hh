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
 * \brief Interface for components that have a liquid state
 */
#ifndef DUMUX_COMPONENT_LIQUID_HH
#define DUMUX_COMPONENT_LIQUID_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief Interface for components that have a liquid state
 */
template<class Scalar, class Component>
class Liquid
{
public:
    /*!
     * \brief Returns true if the liquid phase is assumed to be compressible
     */
    template<class C = Component>
    static constexpr bool liquidIsCompressible()
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: liquidIsCompressible()");
        return 0; // iso c++ requires a return statement for constexpr functions
    }

    /*!
     * \brief Returns true if the liquid phase viscostiy is constant
     */
    template<class C = Component>
    static constexpr bool liquidViscosityIsConstant()
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: liquidViscosityIsConstant()");
        return 0; // iso c++ requires a return statement for constexpr functions
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the liquid component at a given pressure in
     *          \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: liquidDensity(t,p)");
        DUNE_THROW(Dune::NotImplemented, "liquidDensity(t,p)");
    }

    /*!
     * \brief The molar density \f$\mathrm{[mole/m^3]}\f$ of the liquid component at a given pressure in
     *          \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static Scalar liquidMolarDensity(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: liquidMolarDensity(t,p)");
        DUNE_THROW(Dune::NotImplemented, "Component::liquidMolarDensity(t,p)");
    }

    /*!
     * \brief The dynamic liquid viscosity \f$\mathrm{[Pa*s]}\f$ of the pure component.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: liquidViscosity(t,p)");
        DUNE_THROW(Dune::NotImplemented, "liquidViscosity(t,p)");
    }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of the pure component in liquid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static const Scalar liquidEnthalpy(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: liquidEnthalpy(t,p)");
        DUNE_THROW(Dune::NotImplemented, "liquidEnthalpy(t,p)");
    }

    /*!
     * \brief Specific internal energy \f$\mathrm{[J/kg]}\f$ of pure the pure component in liquid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static const Scalar liquidInternalEnergy(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: liquidInternalEnergy(t,p)");
        DUNE_THROW(Dune::NotImplemented, "liquidInternalEnergy(t,p)");
    }

    /*!
     * \brief Thermal conductivity of the component \f$\mathrm{[W/(m*K)]}\f$ as a liquid.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static Scalar liquidThermalConductivity(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: liquidThermalConductivity(t,p)");
        DUNE_THROW(Dune::NotImplemented, "liquidThermalConductivity(t,p)");
    }

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a liquid.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static Scalar liquidHeatCapacity(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: liquidHeatCapacity(t,p)");
        DUNE_THROW(Dune::NotImplemented, "liquidHeatCapacity(t,p)");
    }
};

} // end namespace Components
} // end namespace Dumux

#endif
