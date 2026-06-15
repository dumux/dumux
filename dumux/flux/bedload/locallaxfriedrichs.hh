// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadFlux
 * \copydoc Dumux::Sediment::localLaxFriedrichs
 *
 */
#ifndef DUMUX_FLUX_LOCALLAXFRIEDRICHS_HH
#define DUMUX_FLUX_LOCALLAXFRIEDRICHS_HH

#include "fluxlimiterlet.hh"

namespace Dumux {
namespace Sediment {

/*!
 * \ingroup BedloadFlux
 * \brief Local Lax-Friedrichs like approach to calculate the bedload flux.
 *
 * Apply a Local Lax-Friedrichs like approach to solve the Riemann problem and get the bedload
 * flux \f$q_b\f$ at the cell face.
 *
 * \f[
 * q_b = 0.5 * (q_{b,l} + q_{b,r}) - 0.5 * \lambda_b * (z_r - z_l)
 * \f]
 *
 * where \f$z\f$ is the bed surface. The subscripts \f$l\f$ and \f$r\f$ denote the left
 * the right cell, respectively. \f$\lambda_b\f$ is the bedload wave speed.
 * The bedload wave speed is the maximum of the left and the right cell
 *
 * \f[
 * \lambda_b = max(\lambda_{b,l}, \lambda_{b,r})
 * \f]
 *
 * The bedload wave speed in a cell is estimated using the derivative of the bedload discharge
 * with respect to the bed surface as proposed by Bautista-Parada (2022) in "Decoupled solution
 * of the sediment transport and 2D shallow water equations using the finite volume method".
 * They propose \f$\lambda_b=\xi\frac{3A_g}{h}*u^3\f$, where \f$h\f$ is the water depth,
 * \f$u\f$ the velocity perpendicular to the cell face and \f$\xi=\frac{1}{1-p}\f$ with the porosity \f$p\f$.
 * However, this approach requires the calculation of an equivalent Grass parameter \f$A_g\f$,
 * for each bedload transport formula. Furthermore, the equivalent Grass parameter depends on the friction law,
 * which complicates the implementation of this approach, especially when there are large numbers of bedload
 * transport formulas and friction laws.
 *
 * Therefore, considering the Grass bedload transport formula \f$q_b=A_g*u^3\f$, we substitute \f$A_g\f$
 * by \f$\frac{q_b}{u^3}\f$ which yields:
 *
 * \f[
 * \lambda_b = \xi \frac{3q_b}{h}
 * \f]
 *
 * \param riemannState Contains all quantities needed to define and solve the Riemann state
 *
 * \return bedload flux \f$[m^3/s/m]\f$
 *
 */
template<class Scalar>
Scalar localLaxFriedrichs(const auto& riemannState)
{
    using std::max;
    using std::abs;

    // calculate the bedload wave speed
    Scalar lambdaBedloadLeft = 1 / (1 - riemannState.porosity) * 3 * riemannState.bedloadDischargeLeft / riemannState.waterDepthLeft;
    Scalar lambdaBedloadRight = 1 / (1 - riemannState.porosity) * 3 * riemannState.bedloadDischargeRight / riemannState.waterDepthRight;

    Scalar lambdaBedload = max(abs(lambdaBedloadLeft), abs(lambdaBedloadRight));

    // apply the Local Lax-Friedrichs approach
    Scalar bedloadDischarge =  0.5 * (riemannState.bedloadDischargeLeft + riemannState.bedloadDischargeRight)
                               - 0.5 * lambdaBedload * (riemannState.bedSurfaceRight - riemannState.bedSurfaceLeft);

    // limit the flux based on the sediment availability. This prevents errosion below the fixed bed level
    const Scalar mobility = Sediment::fluxLimiterLET(riemannState.sedimentMassLeft,
                                                     riemannState.sedimentMassRight,
                                                     riemannState.cellAreaLeft,
                                                     riemannState.cellAreaRight,
                                                     bedloadDischarge,
                                                     riemannState.upperMobilityLimitLeft,
                                                     riemannState.lowerMobilityLimitLeft,
                                                     riemannState.upperMobilityLimitRight,
                                                     riemannState.lowerMobilityLimitRight);
    return bedloadDischarge * mobility;
}
} // end namespace Sediment
} // end namespace Dumux

#endif
