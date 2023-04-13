// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief Stress-Strain relationship according to Hooke's law
 *
 * Hooke's law describes a linear stress-strain relationship. \n
 * \n
 * \f$\sigma=\mathbf{c}\epsilon\f$
 * with
 * \f$\epsilon=\frac{1}{2}( \; \mathrm{grad}( \mathbf{u})+\mathrm{grad}(\mathbf{u})^T)\f$ \n
 *\n
 * Furthermore, homogeneous and isotropic materials are assumed. \n
 * As proportionality constants the first and second Lamé parameters are used.\n
 * The final stress-strain relationship yields the following:\n
 * \n
 * \f$ \sigma= 2 \mu \epsilon+\lambda \mathrm{tr}(\epsilon)\mathbf{I} \f$
 * \n
 * \n
 *
 *
 */
#ifndef DUMUX_DISCRETIZATION_HOOKES_LAW_HH
#define DUMUX_DISCRETIZATION_HOOKES_LAW_HH

#include <dumux/flux/hookeslaw_fwd.hh>
#include <dumux/flux/box/hookeslaw.hh>

#endif
