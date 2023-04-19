// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief Advective fluxes according to Forchheims's law (extends Darcy's law by an nonlinear drag)
 *
 * Darcy’s law is linear in the seepage velocity v. As described in "darcyslaw.hh", this is valid for creeping flow (Re<<1).
 * As v increases further, a nonlinear drag arises. The nonlinear drag arises, because the friction becomes comparable to surface drag due to friction.
 * This additional friction is added as the forchheimer term to the Darcy's law. Resulting in the Forchheimer's law.
 * see e.g. Nield & Bejan: Convection in Porous Media \cite nield2006
 *
 * For multiphase flow, the relative passability \f$ \eta_r\f$ is the "Forchheimer-equivalent" to the relative
 * permeability \f$ k_r\f$.
 * We use the same function for \f$ \eta_r\f$ as for \f$ k_r \f$ (Van-Genuchten, Brooks-Corey, linear), other authors use a simple
 * power law e.g.: \f$\eta_{rw} = S_w^3\f$
 *
 * This leads to the equation in the form:
 * \f[ \mathbf{v_\alpha} + c_F \sqrt{\mathbf{K}} \frac{\rho_\alpha}{\mu_\alpha }
 *     |\mathbf{v_\alpha}| \mathbf{v_\alpha}
 *     + \frac{k_{r \alpha}}{\mu_\alpha} \mathbf{K} \nabla \left(p_\alpha
 *     + \rho_\alpha g z \right)=  0
 * \f]
 * This already includes the assumption \f$ k_r(S_w) = \eta_r(S_w)\f$:
 * - \f$\eta_{rw} = S_w^x\f$ looks very similar to e.g. Van Genuchten relative permeabilities
 * - Fichot et al. (2006) \cite Fichot2006 state that several authors claim
 *   that \f$ k_r(S_w), \eta_r(S_w)\f$ can be chosen equal
 * - It leads to the equation not degenerating for the case of \f$S_w=1\f$, because I do not
 *   need to multiply with two different functions, and therefore there are terms not being
 *   zero.
 * - If this assumption is not to be made: Some regularization needs to be introduced ensuring
 *   that not all terms become zero for \f$S_w=1\f$.
 *
 * This non-linear equations is solved for \f$\mathbf{v_\alpha}\f$ using Newton's method
 * and an analytical derivative w.r.t. \f$\mathbf{v_\alpha}\f$.
 *
 * The gradient of the Forchheimer relations looks as follows (mind that \f$\sqrt{\mathbf{K}}\f$
 * is a tensor):
 *
 * \f[  f\left(\mathbf{v_\alpha}\right) =
 * \left(
 * \begin{array}{ccc}
 * 1 & 0 &0 \\
 * 0 & 1 &0 \\
 * 0 & 0 &1 \\
 * \end{array}
 * \right)
 * +
 * c_F \frac{\rho_\alpha}{\mu_\alpha} |\mathbf{v}_\alpha| \sqrt{\mathbf{K}}
 * +
 * c_F \frac{\rho_\alpha}{\mu_\alpha}\frac{1}{|\mathbf{v}_\alpha|} \sqrt{\mathbf{K}}
 * \left(
 * \begin{array}{ccc}
 * v_x^2 & v_xv_y & v_xv_z \\
 * v_yv_x & v_{y}^2 & v_yv_z \\
 * v_zv_x & v_zv_y &v_{z}^2 \\
 * \end{array}
 * \right)
 *  \f]
 *
 * \note We restrict the use of Forchheimer's law to diagonal permeability tensors so far. This might be changed to
 * general tensors using eigenvalue decomposition to get \f$\sqrt{\mathbf{K}}\f$
 * \note Forchheimer's law specialized for different discretization schemes (e.g. Box, CCTpfa).
 * \note Forchheimer's law works with every model that contains Darcy's law.
 */
#ifndef DUMUX_FLUX_FORCHHEIMERS_LAW_HH
#define DUMUX_FLUX_FORCHHEIMERS_LAW_HH

#include <dumux/flux/forchheimerslaw_fwd.hh>

#include <dumux/flux/cctpfa/forchheimerslaw.hh>
#include <dumux/flux/box/forchheimerslaw.hh>

#endif
