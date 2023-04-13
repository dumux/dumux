// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
