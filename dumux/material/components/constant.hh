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
 * \brief Setting constant fluid properties via the input file.
 */
#ifndef DUMUX_COMPONENTS_CONSTANT_HH
#define DUMUX_COMPONENTS_CONSTANT_HH

#include <dune/common/deprecated.hh>
#include <dumux/common/parameters.hh>
#include "component.hh"

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
class Constant : public Component<Scalar, Constant<id, Scalar> >
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
        static const Scalar molarMass = getParamFromGroup<Scalar>(std::to_string(id), "Component.MolarMass", 1.0);
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
        static const Scalar density = getParamFromGroup<Scalar>(std::to_string(id), "Component.LiquidDensity", 1.0);
        return density;
    }

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
        static const Scalar kinematicViscosity = getParamFromGroup<Scalar>(std::to_string(id), "Component.LiquidKinematicViscosity", 1.0);
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
        static const Scalar density = getParamFromGroup<Scalar>(std::to_string(id), "Component.GasDensity", 1.0);
        return density;
    }

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
        static const Scalar kinematicViscosity = getParamFromGroup<Scalar>(std::to_string(id), "Component.GasKinematicViscosity", 1.0);
        return kinematicViscosity * gasDensity(temperature, pressure);
    }
};

} // end namespace Components

} // end namespace Dumux

#endif // DUMUX_COMPONENTS_CONSTANT_HH
