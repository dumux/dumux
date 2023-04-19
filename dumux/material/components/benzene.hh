// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Components
 * \brief A simple benzene component (LNAPL).
 */
#ifndef DUMUX_BENZENE_HH
#define DUMUX_BENZENE_HH

#include <dumux/material/idealgas.hh>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/gas.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A simple benzene component (LNAPL).
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class Benzene
: public Components::Base<Scalar, Benzene<Scalar> >
, public Components::Liquid<Scalar, Benzene<Scalar> >
, public Components::Gas<Scalar, Benzene<Scalar> >
{
     using IdealGas = Dumux::IdealGas<Scalar>;
public:
    /*!
     * \brief A human readable name for the benzene
     */
    static std::string name()
    { return "benzene"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of benzene
     */
    static constexpr Scalar molarMass()
    { return 0.07811; }

    /*!
     * \brief Returns true if the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true if the liquid phase is assumed to be compressible
     */
    static constexpr bool liquidIsCompressible()
    { return false; }

    /*!
     * \brief Returns true if the gas phase viscosity is constant
     */
    static constexpr bool gasViscosityIsConstant()
    { return true; }

    /*!
     * \brief Returns true if the liquid phase viscosity is constant
     */
    static constexpr bool liquidViscosityIsConstant()
    { return true; }

    /*!
     * \brief The density of benzene steam at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
    */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        return IdealGas::density(molarMass(),
                                 temperature,
                                 pressure);
    }

    /*!
     * \brief The molar density of steam in \f$\mathrm{[mol/m^3]}\f$,
     *   depending on pressure and temperature.
     * \param temperature The temperature of the gas
     * \param pressure The pressure of the gas
     */
    static Scalar gasMolarDensity(Scalar temperature, Scalar pressure)
    { return IdealGas::molarDensity(temperature, pressure); }

    /*!
     * \brief The density of pure benzene at a given pressure and temperature \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDensity(Scalar temperature, Scalar pressure)
    {
        return 889.51;
    }

    /*!
     * \brief The molar density of pure benzene at a given pressure and temperature \f$\mathrm{[mol/m^3]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidMolarDensity(Scalar temperature, Scalar pressure)
    {
        return liquidDensity(temperature, pressure)/molarMass();
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure benzene.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidViscosity(Scalar temperature, Scalar pressure)
    {
        return 1.12e-3;//[Pa s]
    }
};

} // end namespace Components

} // end namespace Dumux

#endif
