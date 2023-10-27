// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief Diffusive mass fluxes according to Maxwell-Stefan's law
 *
 * Maxwell-Stefan's law describes the diffusive mass fluxes due to molecular diffusion. The diffusion phenomena results from coupling effects
 * between the different molecules in a gas-mixture \cite Krishna1997. \n
 * The Maxwell-Stefan formulation can be used to describe systems where Fick's law does not hold (e.g. diffusion of diluted
 * gases in multicomponent systems).
 *
 * For diffusive mass fluxes \f$\textbf{j}_{diff}^i\f$ the Maxwell-Stefan formulation can be defined as:
 *
 * \f[
 * \frac{x^i \nabla_T \eta^i}{RT} = - \sum\limits_{j=1,j\neq i}^{N} \frac{x^ix^j}{D^{ij}}\left(\frac{\textbf{j}_{diff}^i}{\varrho^i}-\frac{\textbf{j}_{diff}^j}{\varrho^j}\right) = -
 * \sum\limits_{j=1,j\neq i}^{N} \frac{x^ix^j}{D^{ij}\varrho}\left(\frac{\textbf{j}_{diff}^i}{X^i}-\frac{\textbf{j}_{diff}^j}{X^j}\right)
 * \f]
 *
 * With \f$\eta^i\f$ as the chemical potential of the species i. Note, the diffusion coefficients are based on the Onsager symmetry, thus the diffusion coefficients can be expressed as
 * \f$D^{ij}=D^{ji}\f$.
 *
 */
#ifndef DUMUX_FLUX_MAXWELL_STEFAN_LAW_HH
#define DUMUX_FLUX_MAXWELL_STEFAN_LAW_HH

#include <dumux/flux/maxwellstefanslaw_fwd.hh>

#include <dumux/flux/cctpfa/maxwellstefanslaw.hh>
#include <dumux/flux/box/maxwellstefanslaw.hh>
#include <dumux/flux/staggered/freeflow/maxwellstefanslaw.hh>

#endif
