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
 * \brief  Diffusive heat flux according to non-equilibrium Fourier's law
 *
 * This law is based on the general form of Fourier's law which describes the diffusive
 * heat flux as proportional to a temperature gradient \f$\textbf{grad}\, T_\alpha \f$.
 * In contrast to the general form, a local thermodynamic equilibrium is not assumed.
 * Thus, the heat flux for the different phases \f$\alpha \f$ needs to be solved.
 * \n
 * \f[
 * \textbf{j}_{heat,\alpha} = - \lambda_\alpha \; \textbf{grad}\, T_\alpha
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
