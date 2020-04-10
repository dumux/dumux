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
 * \brief Base class for all components
 * Components provide the thermodynamic relations for the liquid,
 * gaseous and/or solid state of a single
 * chemical species or a _fixed_ mixture of species.
 * Fluid systems use components to compute thermodynamic quantities of phases.
 */
#ifndef DUMUX_COMPONENT_BASE_HH
#define DUMUX_COMPONENT_BASE_HH

#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief IsAqueous struct
 */
template <class Component>
struct IsAqueous : public std::false_type {};

/*!
 * \ingroup Components
 * \brief Base class for all components
 * Components provide the thermodynamic relations for the liquid,
 * gaseous and/or solid state of a single
 * chemical species or a _fixed_ mixture of species.
 * Fluid systems use components to compute thermodynamic quantities of phases.
 */
template <class ScalarType, class Component>
class Base
{
public:

    //! export the scalar type used by the component
    using Scalar = ScalarType;

    //! if the component relies on tabulated values
    static constexpr bool isTabulated = false;

    /*!
     * \brief A default routine for initialization, not needed for components and must not be called.
     *
     * \param tempMin The minimum of the temperature range in \f$\mathrm{[K]}\f$
     * \param tempMax The maximum of the temperature range in \f$\mathrm{[K]}\f$
     * \param nTemp The number of entries/steps within the temperature range
     * \param pressMin The minimum of the pressure range in \f$\mathrm{[Pa]}\f$
     * \param pressMax The maximum of the pressure range in \f$\mathrm{[Pa]}\f$
     * \param nPress The number of entries/steps within the pressure range
     *
     * This function throws a warning when called: "No init routine defined - make sure that this is not necessary!"
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {   Dune::dwarn << "No init routine defined - make sure that this is not necessary!" << std::endl; }

    /*!
     * \brief A human readable name for the component.
     * \note Mandatory for all components
     */
    template<class C = Component>
    static std::string name()
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: name()");
        DUNE_THROW(Dune::NotImplemented, "name()");
    }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of the component.
     */
    template<class C = Component>
    static constexpr Scalar molarMass()
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: molarMass()");
        return 0; // iso c++ requires a return statement for constexpr functions
    }

    /*!
     * \brief Returns the critical temperature in \f$\mathrm{[K]}\f$ of the component.
     */
    template<class C = Component>
    static constexpr Scalar criticalTemperature()
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: criticalTemperature()");
        return 0; // iso c++ requires a return statement for constexpr functions
    }

    /*!
     * \brief Returns the critical pressure in \f$\mathrm{[Pa]}\f$ of the component.
     */
    template<class C = Component>
    static constexpr Scalar criticalPressure()
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: criticalPressure()");
        return 0; // iso c++ requires a return statement for constexpr functions
    }

    /*!
     * \brief Returns the temperature in \f$\mathrm{[K]}\f$ at the component's triple point.
     */
    template<class C = Component>
    static constexpr Scalar tripleTemperature()
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: tripleTemperature()");
        return 0; // iso c++ requires a return statement for constexpr functions
    }

    /*!
     * \brief Returns the pressure in \f$\mathrm{[Pa]}\f$ at the component's triple point.
     */
    template<class C = Component>
    static constexpr Scalar triplePressure()
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: triplePressure()");
        return 0; // iso c++ requires a return statement for constexpr functions
    }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of the component at a given
     *        temperature in \f$\mathrm{[K]}\f$.
     *
     * \param t temperature of the component in \f$\mathrm{[K]}\f$
     */
    template<class C = Component>
    static Scalar vaporPressure(Scalar t)
    {
        static_assert(AlwaysFalse<C>::value, "Mandatory function not implemented: vaporPressure(t)");
        DUNE_THROW(Dune::NotImplemented, "vaporPressure(t)");
    }

};

} // end namespace Components
} // end namespace Dumux

#endif
