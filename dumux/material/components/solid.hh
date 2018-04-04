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
 * \brief Interface for components that have a solid state
 */
#ifndef DUMUX_COMPONENT_SOLID_HH
#define DUMUX_COMPONENT_SOLID_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {
namespace Components {

template<class Scalar, class Component>
class Solid
{
public:
    /*!
     * \brief the component has a solid state if it derives from Solid
     */
    static constexpr bool hasSolidState()
    { return true; }

    /*!
     * \brief Returns true if the solid phase is assumed to be compressible
     */
    template<class C = Component>
    static constexpr bool solidIsCompressible()
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: solidIsCompressible()");
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
    static Scalar solidDensity(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: solidDensity(t,p)");
        DUNE_THROW(Dune::NotImplemented, "solidDensity(t,p)");
    }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of the pure component in solid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static const Scalar solidEnthalpy(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: solidEnthalpy(t,p)");
        DUNE_THROW(Dune::NotImplemented, "solidEnthalpy(t,p)");
    }

    /*!
     * \brief Specific internal energy \f$\mathrm{[J/kg]}\f$ of the pure component in solid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static const Scalar solidInternalEnergy(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: solidInternalEnergy(t,p)");
        DUNE_THROW(Dune::NotImplemented, "solidInternalEnergy(t,p)");
    }

    /*!
     * \brief Thermal conductivity of the component \f$\mathrm{[W/(m*K)]}\f$ as a solid.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static Scalar solidThermalConductivity(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: solidThermalConductivity(t,p)");
        DUNE_THROW(Dune::NotImplemented, "solidThermalConductivity(t,p)");
    }

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a solid.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    template<class C = Component>
    static Scalar solidHeatCapacity(Scalar temperature, Scalar pressure)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: solidHeatCapacity(t,p)");
        DUNE_THROW(Dune::NotImplemented, "solidHeatCapacity(t,p)");
    }

};

} // end namespace Components
} // end namespace Dumux

#endif
