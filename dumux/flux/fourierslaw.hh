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
 * \brief Diffusive heat flux according to Fourier's law
 *
 * According to Fourier's law, the heat flux is proportional to the gradient of the temperature. \n
 * The proportonality constant is the thermal conductivity \f$\lambda \f$. \n
 * The flux is calculated as:\n
 * \n
 * \f[
 * \textbf{j}_{heat} = - \lambda \; \textbf{grad}\, T
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
