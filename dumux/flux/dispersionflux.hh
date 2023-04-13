// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
