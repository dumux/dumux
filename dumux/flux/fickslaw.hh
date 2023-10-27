// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief Diffusive mass flux according to Fick's law
 *
 *
 * Fick's law describes the diffusive flux of mass as proportional to it's concentration gradient in a given phase, caused by the Brownian molecular motion. \n
 * For a single phase system, the proportionality constant is the molecular diffusion coefficient \f$ D_m \f$.
 *
 * \n
 * \f[
 * \mathbf{j}_{d} = - \varrho D_m \nabla  X
 * \f]
 * \n
 *
 * Extending this to multi-phase, multi-component systems, Fick's law can be expressed as follows:
 * \n
 * \f[
 * \mathbf{j}_{d,\alpha}^\kappa = - \varrho_\alpha D_\alpha^\kappa \nabla  X_\alpha^\kappa
 * \f]
 * \n
 *
 * Here \f$D_\alpha^\kappa\f$ is the molecular diffusion coefficient of component \f$\kappa\f$ in phase \f$\alpha\f$.
 * \n
 * In a porous medium, the actual path lines are tortuous due to the impact of the solid matrix. The tortuosity and the impact of
 * the presence of multiple phases is accounted by using an effective diffusion coefficient \f$D_{pm,\alpha}^\kappa\f$. \n
 * The effective diffusion coefficient is then a function of tortuosity \f$\tau\f$, porosity \f$\phi\f$, saturation \f$S\f$ and the molecular diffusion coefficient \f$D_{m}\f$
 * (\f$D_{pm,\alpha}^\kappa=f(\tau,\phi,S_\alpha,D_m)\f$). \n
 * Models to describe those effects are for example Millington-Quirk \cite millington1961 or Constant-Tortuosity \cite carman1937, \cite bear1972.
 * \n
 */
#ifndef DUMUX_FLUX_FICKS_LAW_HH
#define DUMUX_FLUX_FICKS_LAW_HH

#include <dumux/flux/fickslaw_fwd.hh>

#include <dumux/flux/cctpfa/fickslaw.hh>
#include <dumux/flux/ccmpfa/fickslaw.hh>
#include <dumux/flux/box/fickslaw.hh>
#include <dumux/flux/staggered/freeflow/fickslaw.hh>

#endif
