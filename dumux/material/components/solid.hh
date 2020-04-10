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
 * \brief Interface for components that have a solid state
 */
#ifndef DUMUX_COMPONENT_SOLID_HH
#define DUMUX_COMPONENT_SOLID_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief Interface for components that have a solid state
 */
template<class Scalar, class Component>
class Solid
{
public:
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
     */
    template<class C = Component>
    static Scalar solidDensity(Scalar temperature)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: solidDensity(t)");
        DUNE_THROW(Dune::NotImplemented, "solidDensity(t)");
    }

    /*!
     * \brief Thermal conductivity of the component \f$\mathrm{[W/(m*K)]}\f$ as a solid.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    template<class C = Component>
    static Scalar solidThermalConductivity(Scalar temperature)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: solidThermalConductivity(t)");
        DUNE_THROW(Dune::NotImplemented, "solidThermalConductivity(t)");
    }

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a solid.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    template<class C = Component>
    static Scalar solidHeatCapacity(Scalar temperature)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: solidHeatCapacity(t)");
        DUNE_THROW(Dune::NotImplemented, "solidHeatCapacity(t)");
    }

};

} // end namespace Components
} // end namespace Dumux

#endif
