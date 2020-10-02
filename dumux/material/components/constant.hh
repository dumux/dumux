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
 * \brief Setting constant fluid properties via the input file.
 */
#ifndef DUMUX_COMPONENTS_CONSTANT_HH
#define DUMUX_COMPONENTS_CONSTANT_HH

#include <dumux/common/parameters.hh>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/gas.hh>
#include <dumux/material/components/solid.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A component which returns run time specified values
 *        for all fluid properties.
 *
 * \tparam id  The id used to read from the input file / parametertree
 * \tparam Scalar  The type used for scalar values
 *
 * \note For the constant component with id=1 you would specify the parameters in the input file as follows
 *       \code{.ini}
 *       [1.Component]
 *       MolarMass = 0.018 # kg/mol
 *       \endcode
 * \note If you only have one component you can also omit the "1.".
 */
template<int id, class Scalar>
class Constant
: public Components::Base<Scalar, Constant<id, Scalar> >
, public Components::Liquid<Scalar, Constant<id, Scalar> >
, public Components::Gas<Scalar, Constant<id, Scalar> >
, public Components::Solid<Scalar, Constant<id, Scalar> >
{

public:
    /*!
     * \brief Returns true if the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
    { return false; }

    /*!
     * \brief Returns true if the gas phase viscosity is constant
     */
    static constexpr bool gasViscosityIsConstant()
    { return true; }

    /*!
     * \brief Returns true if the gas phase is assumed to be ideal
     */
    static constexpr bool gasIsIdeal()
    { return true; }

    /*!
     * \brief Returns true if the liquid phase is assumed to be compressible
     */
    static constexpr bool liquidIsCompressible()
    { return false; }

    /*!
     * \brief Returns true if the liquid phase viscosity is constant
     */
    static constexpr bool liquidViscosityIsConstant()
    { return true; }

    /*!
     * \brief A human readable name for the component.
     */
    static const std::string& name()
    {
        static const std::string name = getParamFromGroup<std::string>(std::to_string(id), "Component.Name", "component");
        return name;
    }

    /*!
     * \brief The mass in \f$\mathrm{[kg]}\f$ of one mole of the component.
     */
    static Scalar molarMass()
    {
        static const Scalar molarMass = getParamFromGroup<Scalar>(std::to_string(id), "Component.MolarMass");
        return molarMass;
    }

    /*!
     * \brief Sets the liquid density in \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature phase temperature in \f$\mathrm{[K]}\f$
     * \param pressure phase pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        static const Scalar density = getParamFromGroup<Scalar>(std::to_string(id), "Component.LiquidDensity");
        return density;
    }

    /*!
     * \brief The molar density in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar liquidMolarDensity(Scalar temperature, Scalar pressure)
    { return liquidDensity(temperature, pressure)/molarMass(); }

    /*!
     * \brief Sets the liquid dynamic viscosity in \f$\mathrm{[Pa*s]}\f$.
     *
     * Although the dynamic viscosity \f$\mathrm{[Pa*s]}\f$ is returned,
     * the kinematic viscosity \f$\mathrm{[m^2/s]}\f$ is requested from run time input.
     *
     * \param temperature phase temperature in \f$\mathrm{[K]}\f$
     * \param pressure phase pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        static const Scalar kinematicViscosity = getParamFromGroup<Scalar>(std::to_string(id), "Component.LiquidKinematicViscosity");
        return kinematicViscosity * liquidDensity(temperature, pressure);
    }

    /*!
     * \brief Sets the gas density in \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature phase temperature in \f$\mathrm{[K]}\f$
     * \param pressure phase pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        static const Scalar density = getParamFromGroup<Scalar>(std::to_string(id), "Component.GasDensity");
        return density;
    }

    /*!
     * \brief The molar density in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar gasMolarDensity(Scalar temperature, Scalar pressure)
    { return gasDensity(temperature, pressure)/molarMass(); }


    /*!
     * \brief Sets the gas dynamic viscosity in \f$\mathrm{[Pa*s]}\f$.
     *
     * Although the dynamic viscosity \f$\mathrm{[Pa*s]}\f$ is returned,
     * the kinematic viscosity \f$\mathrm{[m^2/s]}\f$ is requested from run time input.
     *
     * \param temperature phase temperature in \f$\mathrm{[K]}\f$
     * \param pressure phase pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        static const Scalar kinematicViscosity = getParamFromGroup<Scalar>(std::to_string(id), "Component.GasKinematicViscosity");
        return kinematicViscosity * gasDensity(temperature, pressure);
    }

        /*!
     * \brief The density in \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure in
     *          \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidDensity(Scalar temperature)
    {
        static const Scalar density = getParamFromGroup<Scalar>(std::to_string(id), "Component.SolidDensity");
        return density;
    }

    /*!
     * \brief Thermal conductivity of the component \f$\mathrm{[W/(m*K)]}\f$ as a solid.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidThermalConductivity(Scalar temperature)
    {
        static const Scalar solidThermalConductivity = getParamFromGroup<Scalar>(std::to_string(id), "Component.SolidThermalConductivity");
        return solidThermalConductivity;
    }

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a solid.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidHeatCapacity(Scalar temperature)
    {
        static const Scalar solidHeatCapacity = getParamFromGroup<Scalar>(std::to_string(id), "Component.SolidHeatCapacity");
        return solidHeatCapacity;
    }
};

} // end namespace Components

} // end namespace Dumux

#endif
