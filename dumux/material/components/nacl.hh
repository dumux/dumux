// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief Material properties of pure salt \f$NaCl\f$.
 */
#ifndef DUMUX_NACL_HH
#define DUMUX_NACL_HH

#include <dumux/common/exceptions.hh>

#include <cmath>
#include <iostream>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/solid.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the NaCl properties
 */
template <class Scalar>
class NaCl
: public Components::Base<Scalar, NaCl<Scalar> >
, public Components::Solid<Scalar, NaCl<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the NaCl.
     */
    static std::string name()
    {
        return "NaCl";
    }

    /*!
     * \brief The molar mass of NaCl in \f$\mathrm{[kg/mol]}\f$.
     */
    static constexpr Scalar molarMass()
    {
        return 58.4428e-3 ;
    }

    /*!
     * \brief The mass density \f$\mathrm{[kg/m^3]}\f$ of NaCl.
     */
    static Scalar solidDensity(Scalar temperature)
    {
        return 2165.0;
    }

    /*!
     * \brief The mass density \f$\mathrm{[kg/m^3]}\f$ of NaCl.
     */
    static Scalar solidMolarDensity(Scalar temperature)
    {
        return solidDensity(temperature)/molarMass();
    }

    /*!
     * \brief The specific heat capacity \f$\mathrm{[J/kg K]}\f$ of NaCl.
     */
    static Scalar solidHeatCapacity(Scalar temperature)
    {
        return 50.50/molarMass();
    }

    /*!
     * \brief Thermal conductivity of the component \f$\mathrm{[W/(m*K)]}\f$ as a solid.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidThermalConductivity(Scalar temperature)
    {
        return 6.49;
    }
};

} // end namespace Components

} // end namespace Dumux

#endif
