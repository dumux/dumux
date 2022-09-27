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
