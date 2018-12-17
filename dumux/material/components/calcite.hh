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
 * \brief A class for the CaCO3 mineral phase properties
 */
#ifndef DUMUX_CALCITE_HH
#define DUMUX_CALCITE_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/solid.hh>

#include <dumux/material/components/calciumion.hh>
#include <dumux/material/components/carbonateion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the CaCO3 mineral phase properties
 */
template <class Scalar>
class Calcite
: public Components::Base<Scalar, Calcite<Scalar> >
, public Components::Solid<Scalar, Calcite<Scalar> >
{

public:
    using CalciumIon = Components::CalciumIon<Scalar>;
    using CarbonateIon = Components::CarbonateIon<Scalar>;
    /*!
     * \brief A human readable name for calcite.
     */
    static std::string name()
    { return "Calcite"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of calcite.
     */
    static constexpr Scalar molarMass()
    { return CalciumIon::molarMass() + CarbonateIon::molarMass(); } // kg/mol

    /*!
     * \brief Returns true if the solid phase is assumed to be compressible
     */
    static constexpr bool solidIsCompressible()
    { return false; }

    /*!
     * \brief The density in \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure in
     *          \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static constexpr Scalar solidDensity(Scalar temperature)
    { return 2.71e3; }

    /*!
     * \brief Thermal conductivity of the component \f$\mathrm{[W/(m*K)]}\f$ as a solid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidThermalConductivity(Scalar temperature)
    { return 3.849; }

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a solid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidHeatCapacity(Scalar temperature)
    { return 837; }
};

} // end namespace Components
} // end namespace Dumux

#endif
