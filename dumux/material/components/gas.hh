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
 * \brief Interface for components that have a gas state
 */
#ifndef DUMUX_COMPONENT_GAS_HH
#define DUMUX_COMPONENT_GAS_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief Interface for components that have a gas state
 */
template<class Scalar, class Component>
class Gas
{
public:
    /*!
     * \brief Returns true if the gas phase is assumed to be compressible
     */
    template<class C = Component>
    static constexpr bool gasIsCompressible()
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: gasIsCompressible()");
        return 0; // iso c++ requires a return statement for constexpr functions
    }

    /*!
     * \brief Returns true if the gas phase viscostiy is constant
     */
    template<class C = Component>
    static constexpr bool gasViscosityIsConstant()
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: gasViscosityIsConstant()");
        return 0; // iso c++ requires a return statement for constexpr functions
    }

    /*!
     * \brief Returns true if the gas phase is assumed to be ideal
     */
    template<class C = Component>
    static constexpr bool gasIsIdeal()
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: gasIsIdeal()");
        return 0; // iso c++ requires a return statement for constexpr functions
    }

    /*!
     * \brief The density in \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure in
     *          \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: gasDensity(t,p)");
        DUNE_THROW(Dune::NotImplemented, "gasDensity(t,p)");
    }

    /*!
     * \brief The molar density in \f$\mathrm{[mol/m^3]}\f$ of the component at a given pressure in
     *          \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static Scalar gasMolarDensity(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: gasMolarDensity(t,p)");
        DUNE_THROW(Dune::NotImplemented, "gasMolarDensity(t,p)");
    }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of the pure component in gas.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static const Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: gasEnthalpy(t,p)");
        DUNE_THROW(Dune::NotImplemented, "gasEnthalpy(t,p)");
    }

    /*!
     * \brief Specific internal energy \f$\mathrm{[J/kg]}\f$ of the pure component in gas.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static const Scalar gasInternalEnergy(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: gasInternalEnergy(t,p)");
        DUNE_THROW(Dune::NotImplemented, "gasInternalEnergy(t,p)");
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of the pure component at a given pressure in
     *          \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: gasViscosity(t,p)");
        DUNE_THROW(Dune::NotImplemented, "gasViscosity(t,p)");
    }

    /*!
     * \brief Thermal conductivity of the component \f$\mathrm{[W/(m*K)]}\f$ as a gas.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static Scalar gasThermalConductivity(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: gasThermalConductivity(t,p)");
        DUNE_THROW(Dune::NotImplemented, "gasThermalConductivity(t,p)");
    }

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a gas.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static Scalar gasHeatCapacity(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: gasHeatCapacity(t,p)");
        DUNE_THROW(Dune::NotImplemented, "gasHeatCapacity(t,p)");
    }

};

} // end namespace Components
} // end namespace Dumux

#endif
