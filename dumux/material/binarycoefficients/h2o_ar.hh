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
        static const Scalar E = 2310.5463;
        static const Scalar F = -46.7034;
        static const Scalar G = 160.4066;
        static const Scalar H = -118.3043;

        return henryIAPWS(E, F, G, H, temperature);
    }

};

} // end namespace BinaryCoeff
} // end namespace Dumux

#endif
