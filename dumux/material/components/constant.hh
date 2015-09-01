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
 * \brief Setting constant fluid properties via the input file for testing purposes.
 */
#ifndef DUMUX_CONSTANT_HH
#define DUMUX_CONSTANT_HH

#include <dumux/common/parameters.hh>

#include "component.hh"

namespace Dumux
{
/*!
 * \ingroup Components
 *
 * \brief A component which returns run time specified values
 *        for all fluid properties.
 *
 * \tparam Scalar  The type used for scalar values
 */
template<class TypeTag, class Scalar>
class Constant : public Component<Scalar, Constant<TypeTag, Scalar> >
{

public:
    /*!
     * \brief A human readable name for the component.
     */
    static const char *name()
    { return "Constant"; }

    /*!
     * \brief Returns true if the liquid phase is assumed to be compressible
     */
    static bool liquidIsCompressible()
    { return false; }

    /*!
     * \brief Sets the liquid density in \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature phase temperature in \f$\mathrm{[K]}\f$
     * \param pressure phase pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        static const Scalar density
            = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, LiquidDensity);
        return density;
    }

    /*!
     * \brief Sets the liquid dynamic viscosity in \f$\mathrm{[Pa*s]}\f$.
     *
     * Although the dynamic viscosity \f$\mathrm{[Pa*s]} is returned,
     * the kinematic viscosity \f$\mathrm{[m^2/s]} is requested from run time input.
     *
     * \param temperature phase temperature in \f$\mathrm{[K]}\f$
     * \param pressure phase pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        static const Scalar kinematicViscosity
            = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, LiquidKinematicViscosity);
        return kinematicViscosity * liquidDensity(temperature, pressure);
    }


    /*!
     * \brief Returns true if the gas phase is assumed to be compressible
     */
    static bool gasIsCompressible()
    { return false; }

    /*!
     * \brief Sets the gas density in \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature phase temperature in \f$\mathrm{[K]}\f$
     * \param pressure phase pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        static const Scalar density
            = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, GasDensity);
        return density;
    }

    /*!
     * \brief Sets the gas dynamic viscosity in \f$\mathrm{[Pa*s]}\f$.
     *
     * Although the dynamic viscosity \f$\mathrm{[Pa*s]} is returned,
     * the kinematic viscosity \f$\mathrm{[m^2/s]} is requested from run time input.
     *
     * \param temperature phase temperature in \f$\mathrm{[K]}\f$
     * \param pressure phase pressure in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        static const Scalar kinematicViscosity
            = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, GasKinematicViscosity);
        return kinematicViscosity * gasDensity(temperature, pressure);
    }
};

} // end namespace

#endif // DUMUX_CONSTANT_HH
