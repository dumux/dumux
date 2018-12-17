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
 * \brief Properties of pure molecular oxygen \f$O_2\f$.
 */
#ifndef DUMUX_GRANITE_HH
#define DUMUX_GRANITE_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/solid.hh>

#include <cmath>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief Properties of granite
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class Granite : public Components::Base<Scalar, Granite<Scalar> >
              , public Components::Solid<Scalar, Granite<Scalar> >

{
public:
    /*!
     * \brief A human readable name for the solid.
     */
    static std::string name()
    { return "Granite"; }

    /*!
     * \brief Returns true if the solid phase is assumed to be compressible
     */
    static constexpr bool solidIsCompressible()
    {
        return false; // iso c++ requires a return statement for constexpr functions
    }

    /*!
     * \brief The molar mass of Siliciumoxide which is 70 % of granite in \f$\mathrm{[kg/mol]}\f$.
     */
    static constexpr Scalar molarMass()
    {
        return 60.08e-3;
    }

    /*!
     * \brief The density in \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure in
     *          \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidDensity(Scalar temperature)
    {
       return 2700;
    }

    /*!
     * \brief Thermal conductivity of the component \f$\mathrm{[W/(m*K)]}\f$ as a solid.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidThermalConductivity(Scalar temperature)
    {
        return 2.8;
    }

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a solid.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidHeatCapacity(Scalar temperature)
    {
        return 790;
    }

};

} // end namespace Components

} // end namespace Dumux

#endif
