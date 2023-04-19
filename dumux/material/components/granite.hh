// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
