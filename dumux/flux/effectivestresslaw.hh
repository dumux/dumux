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
 * \brief Effective stress are used to describe the actual stresses acting on the grains/matrix in the soil.
 * Furthermore, they determine the behaviour of the soil.
 * Most of the geomechanical laws are written in terms of effective stresses,
 * such as the poroelastic model in dumux/geomechanics/poroelastic/model.hh
 *
 * Effective stresses are denoted with \f$ \boldsymbol{\sigma_{\mathrm{eff}}} \f$.
 * The effective stress tensor \f$ \boldsymbol{\sigma_{\mathrm{eff}}} \f$ is
 * determined by the stress tensor \f$ \boldsymbol{\sigma} \f$ , the effective pore pressure \f$ p_{\mathrm{eff}} \f$ and the Biot's coefficient \f$ \alpha \f$ :
 * \f[
 * \boldsymbol{\sigma_{\mathrm{eff}}} = \boldsymbol{\sigma} - \alpha p_{\mathrm{eff}} \mathbf{I}
 * \f]
 * \note The implementations of the effective stress laws are inside the folders of their respective discretization.
 *
 *
 */
#ifndef DUMUX_FLUX_EFFECIVESTRESS_LAW_HH
#define DUMUX_FLUX_EFFECIVESTRESS_LAW_HH

#include <dumux/flux/effectivestresslaw_fwd.hh>
#include <dumux/flux/box/effectivestresslaw.hh>

#endif
