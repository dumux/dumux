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
#ifndef DUMUX_BINARY_COEFF_H2O_AR_HH
#define DUMUX_BINARY_COEFF_H2O_AR_HH

#include <dumux/material/binarycoefficients/henryiapws.hh>

namespace Dumux {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and argon.
 */
class H2O_AR
{
public:
  /*!
     * \brief Henry coefficient \f$\mathrm{[Pa]}\f$  for argon in liquid water.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     */
    template <class Scalar>
    static Scalar henry(Scalar temperature)
    {
        static constexpr Scalar E = 2310.5463;
        static constexpr Scalar F = -46.7034;
        static constexpr Scalar G = 160.4066;
        static constexpr Scalar H = -118.3043;

        return henryIAPWS(E, F, G, H, temperature);
    }

    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure) {
        DUNE_THROW(Dune::NotImplemented, "BinaryCoefficients::H2O_AR::gasDiffCoeff()");
    }

    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure) {
        DUNE_THROW(Dune::NotImplemented, "BinaryCoefficients::H2O_AR::liquidDiffCoeff()");
    }

    template <class Scalar>
    static Scalar henryMixture(Scalar temperature) {
        DUNE_THROW(Dune::NotImplemented, "BinaryCoefficients::H2O_AR::henryMixture()");
    }

};

} // end namespace BinaryCoeff
} // end namespace Dumux

#endif
