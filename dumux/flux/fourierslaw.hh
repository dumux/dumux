// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief Diffusive heat flux according to Fourier's law
 *
 * According to Fourier's law, the heat flux is proportional to the gradient of the temperature. \n
 * The proportonality constant is the thermal conductivity \f$\lambda \f$. \n
 * The flux is calculated as:\n
 * \n
 * \f[
 * \textbf{j}_{heat} = - \lambda \; \nabla  T
 * \f]
 * \n
 * \n
 * In general, for porous medium local thermal equilibrium is assumed.\n
 * To account for effects of porous media, an effective thermal conductivity \f$\lambda_{pm}\f$ is used.\n
 * Models for the effective thermal conductivity are for example: Johansen \cite johansen1977,  simple fluid lumping, Somerton \cite somerton1974.
 *
 */
#ifndef DUMUX_FLUX_FOURIERS_LAW_HH
#define DUMUX_FLUX_FOURIERS_LAW_HH

#include <dumux/flux/fourierslaw_fwd.hh>

#include <dumux/flux/cctpfa/fourierslaw.hh>
#include <dumux/flux/ccmpfa/fourierslaw.hh>
#include <dumux/flux/box/fourierslaw.hh>
#include <dumux/flux/staggered/freeflow/fourierslaw.hh>

#endif
