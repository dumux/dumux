// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Binarycoefficients
 * \brief Various relations for molecular diffusion coefficients.
 */
#ifndef DUMUX_FULLERMETHOD_HH
#define DUMUX_FULLERMETHOD_HH

#include <dumux/common/math.hh>

namespace Dumux {
namespace BinaryCoeff {
/*!
 * \ingroup Binarycoefficients
 * \brief Estimate binary diffusion coefficients \f$\mathrm{[m^2/s]}\f$ in gases according to
 *        the method by Fuller.
 *
 * \param M molar masses \f$\mathrm{[g/mol]}\f$
 * \param SigmaNu atomic diffusion volume
 * \param temperature The temperature \f$\mathrm{[K]}\f$
 * \param pressure phase pressure \f$\mathrm{[Pa]}\f$
 *
 * This function estimates the diffusion coefficients in binary gases
 * using to the method proposed by Fuller. This method and is only
 * valid at "low" pressures.
 *
 * See: R. Reid, et al. (1987, pp. 587-588) \cite reid1987
 */
template <class Scalar>
inline Scalar fullerMethod(const Scalar *M, // molar masses [g/mol]
                           const Scalar *SigmaNu, // atomic diffusion volume
                           const Scalar temperature, // [K]
                           const Scalar pressure) // [Pa]
{
    // "effective" molar mass in [g/m^3]
    Scalar Mab = harmonicMean(M[0], M[1]);

    // Fuller's method
    using std::pow;
    using std::sqrt;
    Scalar tmp = pow(SigmaNu[0], 1./3) + pow(SigmaNu[1], 1./3);
    return 1e-4 * (143.0*pow(temperature, 1.75))/(pressure*sqrt(Mab)*tmp*tmp);
}

} // end namespace BinaryCoeff
} // end namespace Dumux

#endif
