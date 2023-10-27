// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief  Diffusive heat flux according to non-equilibrium Fourier's law
 *
 * This law is based on the general form of Fourier's law which describes the diffusive
 * heat flux as proportional to a temperature gradient \f$\nabla  T_\alpha \f$.
 * In contrast to the general form, a local thermodynamic equilibrium is not assumed.
 * Thus, the heat flux for the different phases \f$\alpha \f$ needs to be solved.
 * \n
 * \f[
 * \textbf{j}_{heat,\alpha} = - \lambda_\alpha \; \nabla  T_\alpha
 * \f]
 * \n
 * With \f$\lambda_\alpha \f$ as the thermal conductivity for either a solid, liquid or
 * gaseous phase.
 */
#ifndef DUMUX_DISCRETIZATION_FOURIERS_LAW_NONEQUILIBRIUM_HH
#define DUMUX_DISCRETIZATION_FOURIERS_LAW_NONEQUILIBRIUM_HH

#include <dumux/flux/fourierslawnonequilibrium_fwd.hh>

#include <dumux/flux/cctpfa/fourierslawnonequilibrium.hh>
#include <dumux/flux/box/fourierslawnonequilibrium.hh>

#endif
