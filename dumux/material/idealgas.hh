// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Material
 * \brief Relations valid for an ideal gas.
 */
#ifndef DUMUX_IDEAL_GAS_HH
#define DUMUX_IDEAL_GAS_HH

#include <dumux/material/constants.hh>

namespace Dumux {

/*!
 * \ingroup Material
 * \brief Relations valid for an ideal gas.
 */
template <class Scalar>
class IdealGas
{
public:
    //! The ideal gas constant \f$\mathrm{[J/mol/K]}\f$
    static constexpr Scalar R = Constants<Scalar>::R;

    /*!
     * \brief The density of the gas in \f$\mathrm{[kg/m^3]}\f$, depending on
     *        pressure, temperature and average molar mass of the gas.
     * \param avgMolarMass The average molar mass of the gas
     * \param temperature The temperature of the gas
     * \param pressure The pressure of the gas
     */
    static constexpr Scalar density(Scalar avgMolarMass,
                                    Scalar temperature,
                                    Scalar pressure)
    { return molarDensity(temperature,pressure)*avgMolarMass;}

    /*!
     * \brief The pressure of the gas in \f$\mathrm{[Pa]}\f$, depending on
     *        the molar density and temperature.
     * \param temperature The temperature of the gas
     * \param rhoMolar The molar density of the gas
     */
    static constexpr Scalar pressure(Scalar temperature,
                                     Scalar rhoMolar)
    { return R*temperature*rhoMolar; }

    /*!
     * \brief The molar density of the gas \f$\mathrm{[mol/m^3]}\f$,
     *        depending on pressure and temperature.
     * \param temperature The temperature of the gas
     * \param pressure The pressure of the gas
     */
    static constexpr Scalar molarDensity(Scalar temperature,
                                         Scalar pressure)
    { return pressure/(R*temperature); }
};
} // end namespace

#endif
