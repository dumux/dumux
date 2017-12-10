// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Calculates the residual of models based on the box scheme element-wise.
 */
#ifndef DUMUX_NAVIERSTOKES_NC_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_NC_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/freeflow/navierstokes/localresidual.hh>
#include <dumux/freeflow/navierstokesnc/staggered/localresidual.hh>

namespace Dumux
{

/*!
 *
 * \todo Please doc me more!
 */

// // forward declaration
template<class TypeTag, DiscretizationMethods Method>
class NavierStokesNCResidualImpl;

template<class TypeTag>
using NavierStokesNCResidual = NavierStokesNCResidualImpl<TypeTag, GET_PROP_VALUE(TypeTag, DiscretizationMethod)>;

}

#endif   // DUMUX_NAVIERSTOKES_NC_LOCAL_RESIDUAL_HH
