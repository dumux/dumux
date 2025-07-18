// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 *
 * \brief A single-phase, isothermal Navier-Stokes model
 *
 * This model implements a single-phase, isothermal Navier-Stokes model, solving the <B> momentum balance equation </B>
 * \f[
 \frac{\partial (\varrho \textbf{v})}{\partial t} + \nabla \cdot (\varrho \textbf{v} \textbf{v}^{\text{T}}) = \nabla \cdot (\mu (\nabla \textbf{v} + \nabla \textbf{v}^{\text{T}}))
   - \nabla p + \varrho \textbf{g} - \textbf{f}
 * \f]
 * By setting the runtime parameter <code>Problem.EnableInertiaTerms</code> to <code>false</code> the Stokes
 * equation can be solved. In this case the term
 * \f[
 *    \nabla \cdot (\varrho \textbf{v} \textbf{v}^{\text{T}})
 * \f]
 * is neglected.
 *
 * The <B> mass balance equation </B>
 * \f[
       \frac{\partial \varrho}{\partial t} + \nabla \cdot (\varrho \textbf{v}) - q = 0
 * \f]
 *
 * closes the system.
 *
 */

#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_MODEL_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_MODEL_HH

#include <dumux/freeflow/navierstokes/momentum/fcstaggered/model.hh>
#warning "This header is deprecated and will be removed after 3.11. Use dumux/freeflow/navierstokes/momentum/fcstaggered/model.hh instead."

#endif // DUMUX_NAVIERSTOKES_MODEL_HH
