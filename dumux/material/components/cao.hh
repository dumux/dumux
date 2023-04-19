// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief Material properties of pure Calcium-Oxide \f$CaO\f$.
 */
#ifndef DUMUX_CAO_HH
#define DUMUX_CAO_HH

#include <dumux/common/exceptions.hh>

#include <cmath>
#include <iostream>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/solid.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the CaO properties
 */
template <class Scalar>
class CaO
: public Components::Base<Scalar, CaO<Scalar> >
, public Components::Solid<Scalar, CaO<Scalar> >
{
public:
    /*!
     * \brief A human readable name for the CaO.
     */
    static const char *name()
    {
        return "CaO";
    }

    /*!
     * \brief The molar mass of CaOH2 in \f$\mathrm{[kg/mol]}\f$.
     */
    static constexpr Scalar molarMass()
    {
        return 56.0774e-3;
    }

    /*!
     * \brief The mass density \f$\mathrm{[kg/m^3]}\f$ of CaO.
     */
    static Scalar solidDensity(Scalar temperature)
    {
        return 3370;
    }

    /*!
     * \brief The molar density \f$\mathrm{[mol/m^3]}\f$ of CaO.
     * Molar density at 293 K. Literature value from Shao et al. (2013).
     */
    static Scalar solidMolarDensity(Scalar temperature)
    {
        return solidDensity(temperature)/molarMass();
    }

    /*!
     * \brief The specific heat capacity \f$\mathrm{[J/kg K]}\f$ of CaO.
     */
    static Scalar solidHeatCapacity(Scalar temperature)
    {
        return 934;  //Nagel et al. (2014) : 934 J/kgK
    }

    /*!
     * \brief The thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     */
    static Scalar solidThermalConductivity(Scalar temperature)
    {
        return  0.4;  //Nagel et al. (2014)
    }
};

} // end namespace Components

} // end namespace Dumux

#endif
