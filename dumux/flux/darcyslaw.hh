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
 * \brief Advective fluxes according to Darcy's law
 *
 * Darcy's law describes the advective flux in porous media on the macro-scale and is valid in the creeping flow regime (Reynolds number << 1).
 * The advective flux characterizes the bulk flow for each fluid phase including all components in case of compositional flow.
 * It is driven by the potential gradient \f$\textbf{grad}\, p - \varrho {\textbf g}\f$,
 * accounting for both pressure-driven and gravitationally-driven flow.
 * The velocity is proportional to the potential gradient with the proportional factor \f$\frac{\textbf K}{\mu}\f$,
 * including the intrinsic permeability of the porous medium, and the viscosity µ of the fluid phase. For one-phase flow it is:
 * \f[
 * v = - \frac{\mathbf K}{\mu}
 * \left(\textbf{grad}\, p - \varrho {\mathbf g} \right)
 * \f]
 * This equation can be extended to calculate the velocity \f$v_\alpha\f$ of phase \f$\alpha\f$ in the case of multi-phase
 * flow by introducing a relative permeability \f$k_{r\alpha}\f$ restricting flow in the presence of other phases:
 * \f[
 * v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 * \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right)
 * \f]
 *
 * Darcy's law is specialized for different discretization schemes.
 * This file contains the data which is required to calculate
 * volume and mass fluxes of fluid phases over a face of a finite volume by means
 * of the Darcy approximation. See the corresponding header files for the specific different discretization methods.
 */
#ifndef DUMUX_FLUX_DARCYS_LAW_HH
#define DUMUX_FLUX_DARCYS_LAW_HH

#include <dumux/flux/darcyslaw_fwd.hh>

#include <dumux/flux/cctpfa/darcyslaw.hh>
#include <dumux/flux/ccmpfa/darcyslaw.hh>
#include <dumux/flux/cvfe/darcyslaw.hh>

#endif
