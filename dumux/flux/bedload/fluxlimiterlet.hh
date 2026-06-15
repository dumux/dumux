// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadFlux
 * \copydoc Dumux::Sediment::fluxLimiterLET
 *
 */
#ifndef DUMUX_FLUX_LIMITER_LET_HH
#define DUMUX_FLUX_LIMITER_LET_HH

namespace Dumux {
namespace Sediment {

/*!
 * \ingroup BedloadFlux
 * \brief Provide a mobility to scale the bedload flux for small sediment masses.
 *
 * This function calculates a mobility, which can be multiplied with the bedload flux to scale it.
 * The mobility depends on the sediment mass in the upwind cell. If the sediment mass is above an upper limit \f$m_{upper}\f$,
 * the mobility is 1 and therefore no limiting applied. If the sediment mass is below a lower limit \f$m_{lower}\f$, the mobility is 0.
 * In this case the bedload discharge is set to zero.
 * Between the upper and the lower limit the mobility is scaled using a modiefied version of the LET model
 * presented by Lomeland et al (2005) in "A New Versatile Relative Permeability Correlation".
 * The mobility is calculated as suggested by Lomeland et al:
 * \f[
 * mobility = \frac{k_{rw} * m_{norm}^L}{m_{norm}^L + E * (1.0 - m_{norm})^T}
 * \f]
 * the parameter \f$m_{norm}\f$ is calculated differently. Since the original approach is made to limit saturation values
 * and can deal therefore only with values smaler than 1.0.
 * \f[
 * m_{norm} = max\left(min\left(\frac{m - m_{lower}}{m_ {upper}-m_{lower}}, 1.0\right), 0.0\right)
 * \f]
 * The parameters are fixed to \f$k_{rw}=1.0\f$, \f$L=2.0\f$, \f$T=1.5\f$ and \f$E=2.0\f$ to ensure a good limiting behaviour.
 * They have no physical meaning.
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \param sedimentMassLeft The sediment mass in the left cell \f$[kg]\f$
 * \param sedimentMassRight The sediment mass in the right cell \f$[kg]\f$
 * \param cellAreaLeft The area of the left cell \f$[m^2]\f$
 * \param cellAreaRight The area of the right cell \f$[m^2]\f$
 * \param localFlux flux at the interface  \f$[m^3/s/m]\f$
 * \param upperMobilityLimitLeft Upper mobility limit of the left cell \f$[kg/m^2]\f$
 * \param lowerMobilityLimitLeft Lower mobility limit of the left cell \f$[kg/m^2]\f$
 * \param upperMobilityLimitRight Upper mobility limit of the right cell \f$[kg/m^2]\f$
 * \param lowerMobilityLimitRight Lower mobility limit of the right cell \f$[kg/m^2]\f$
 *
 * \return The mobility (0 <= mobility <= 1) \f$[-]\f$
 */
template<class Scalar>
static Scalar fluxLimiterLET(const Scalar sedimentMassLeft,
                             const Scalar sedimentMassRight,
                             const Scalar cellAreaLeft,
                             const Scalar cellAreaRight,
                             const Scalar localFlux,
                             const Scalar upperMobilityLimitLeft,
                             const Scalar lowerMobilityLimitLeft,
                             const Scalar upperMobilityLimitRight,
                             const Scalar lowerMobilityLimitRight)
{
    using std::pow;
    using std::clamp;

    Scalar sedimentMass;
    Scalar cellArea;
    Scalar upperMobilityLimit;
    Scalar lowerMobilityLimit;

    // Use an upwind method to determine the cell relevant for the mobility calculation.
    if (localFlux < 0)
    {
        sedimentMass = sedimentMassRight;
        cellArea = cellAreaRight;
        upperMobilityLimit = upperMobilityLimitRight;
        lowerMobilityLimit = lowerMobilityLimitRight;
    }
    else
    {
        sedimentMass = sedimentMassLeft;
        cellArea = cellAreaLeft;
        upperMobilityLimit = upperMobilityLimitLeft;
        lowerMobilityLimit = lowerMobilityLimitLeft;
    }

    // normalize the sediment mass between the mobility limits
    const Scalar m_norm = clamp((sedimentMass - lowerMobilityLimit * cellArea) / (upperMobilityLimit * cellArea - lowerMobilityLimit * cellArea), 0.0, 1.0);

    // use the LET-model to get the mobility
    constexpr Scalar krw = 1.0;  // default value is 1.0
    constexpr Scalar letL = 2.0;  // default value is 2.0
    constexpr Scalar letT = 1.5;  // default value is 2.0
    constexpr Scalar letE = 2.0;  // default value is 1.0
    return krw * pow(m_norm, letL)/(pow(m_norm, letL) + letE * pow(1.0 - m_norm, letT));
}

} // end namespace Sediment
} // end namespace Dumux

#endif
