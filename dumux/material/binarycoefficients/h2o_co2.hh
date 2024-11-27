// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and oxygen.
 */
#ifndef DUMUX_BINARY_COEFF_H2O_CO2_HH
#define DUMUX_BINARY_COEFF_H2O_CO2_HH

#include <dumux/material/binarycoefficients/henryiapws.hh>

namespace Dumux {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and carbon dioxide.
 */
class H2O_CO2
{
public:
  /*!
     * \brief Henry coefficient \f$\mathrm{[Pa]}\f$  for carbon dioxide in liquid water.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     */
    template <class Scalar>
    static Scalar henry(Scalar temperature)
    {
        static constexpr Scalar E = 1672.9376;
        static constexpr Scalar F = 28.1751;
        static constexpr Scalar G = -112.4619;
        static constexpr Scalar H = 85.3807;

        return henryIAPWS(E, F, G, H, temperature);
    }

    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure) {
        DUNE_THROW(Dune::NotImplemented, "BinaryCoefficients::H2O_CO2::gasDiffCoeff()");
    }

    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure) {
        DUNE_THROW(Dune::NotImplemented, "BinaryCoefficients::H2O_CO2::liquidDiffCoeff()");
    }

    template <class Scalar>
    static Scalar henryMixture(Scalar temperature) {
        DUNE_THROW(Dune::NotImplemented, "BinaryCoefficients::H2O_CO2::henryMixture()");
    }

};

} // end namespace BinaryCoeff
} // end namespace Dumux

#endif
