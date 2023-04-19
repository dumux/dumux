// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for nitrogen and oxygen.
 */
#ifndef DUMUX_BINARY_COEFF_N2_O2_HH
#define DUMUX_BINARY_COEFF_N2_O2_HH

#include <dumux/material/binarycoefficients/henryiapws.hh>
#include <dumux/material/binarycoefficients/fullermethod.hh>

#include <dumux/material/components/o2.hh>
#include <dumux/material/components/n2.hh>

namespace Dumux {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for nitrogen and oxygen.
 */
class N2_O2
{
public:
    /*!
     * \brief Henry coefficient \f$\mathrm{[Pa]}\f$  for molecular oxygen in liquid nitrogen.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     */
    template <class Scalar>
    static Scalar henry(Scalar temperature)
    {
        DUNE_THROW(Dune::NotImplemented, "henry coefficient for oxygen in liquid nitrogen");
    }

    /*!
     * \brief Binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular oxygen in liquid nitrogen.
     *
     * Uses fullerMethod to determine the diffusion of water in nitrogen.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        using N2 = Dumux::Components::N2<Scalar>;
        using O2 = Dumux::Components::O2<Scalar>;

        // atomic diffusion volumes
        const Scalar SigmaNu[2] = { 18.1 /* N2 */,  16.3 /* O2 */ };
        // molar masses [g/mol]
        const Scalar M[2] = { N2::molarMass()*1e3, O2::molarMass()*1e3 };
        return fullerMethod(M, SigmaNu, temperature, pressure);
    }

    /*!
     * \brief Diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular oxygen in liquid nitrogen.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "diffusion coefficient for liquid oxygen and nitrogen");
    }
};

} // end namespace BinaryCoeff
} // end namespace Dumux

#endif
