// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and heavy oil.
 */
#ifndef DUMUX_BINARY_COEFF_H2O_HEAVYOIL_HH
#define DUMUX_BINARY_COEFF_H2O_HEAVYOIL_HH

#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/heavyoil.hh>

namespace Dumux {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and heavy oil as in SAGD processes
 */
class H2O_HeavyOil
{
public:
    /*!
     * \brief Henry coefficient \f$[N/m^2]\f$  for heavy oil in liquid water.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     *
     * \todo values copied from TCE, please improve it
     */
    template <class Scalar>
    static Scalar henryOilInWater(Scalar temperature)
    {
        // values copied from TCE, TODO: improve this!!
        Scalar dumuxH = 1.5e-1 / 101.325; // unit [(mol/m^3)/Pa]
        dumuxH *= 18.02e-6;  //multiplied by molar volume of reference phase = water
        return 1.0/dumuxH; // [Pa]
    }

    /*!
     * \brief Henry coefficient \f$[N/m^2]\f$  for water in liquid heavy oil.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     *
     * \todo arbitrary value, please improve it
     */
    template <class Scalar>
    static Scalar henryWaterInOil(Scalar temperature)
    {
        // arbitrary, TODO: improve it!!
        return 1.0e8; // [Pa]
    }


    /*!
     * \brief Binary diffusion coefficient [m^2/s] for molecular water and heavy oil.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     *
     * \todo value is just an order of magnitude, please improve it
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        return 1e-6; // [m^2/s] TODO: This is just an order of magnitude. Please improve it!
    }

    /*!
     * \brief Diffusion coefficient [m^2/s] for heavy oil in liquid water.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     *
     * \todo value is just an order of magnitude, please improve it
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        return 1.e-9;  // This is just an order of magnitude. Please improve it!
    }
};

} // end namespace BinaryCoeff
} // end namespace Dumux

#endif
