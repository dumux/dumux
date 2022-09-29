// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Flux
 * \brief Dispersion flux for different discretization schemes
 *
 * Phenomenologically, dispersion is a similar process to diffusion. However, the additional spreading of components is due to fluctuations of magnitude and direction of the flow
velocities. \n
By replacing the diffusion coefficient \f$D_\alpha^\kappa\f$ with \f$D_{\alpha,eff}^\kappa\f$ in Fick's law the velocity-dependent effects of dispersion can be expressed \cite bear1972.
 *
 * \f[
 * D_{\alpha,eff}^\kappa = D_\alpha^\kappa + D_{\alpha,disp}^\kappa(\textbf{v}_\alpha)
 * \f]
 *
 * Possible options of describing the dispersion tensors can be found in \cite scheidegger1961.
 */
#ifndef DUMUX_FLUX_DISPERSION_FLUX_HH
#define DUMUX_FLUX_DISPERSION_FLUX_HH

#include <dumux/flux/dispersionflux_fwd.hh>
#include <dumux/flux/box/dispersionflux.hh>
#include <dumux/flux/cctpfa/dispersionflux.hh>

#endif
