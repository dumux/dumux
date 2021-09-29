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
 * \ingroup SpatialParameters
 * \brief A spatial params implementation for 1p problem with constant properties
 */
#ifndef DUMUX_FV_CONSTANT_SPATIAL_PARAMS_ONE_P_HH
#define DUMUX_FV_CONSTANT_SPATIAL_PARAMS_ONE_P_HH

#warning "This header will be removed after 3.5 in favor of dumux/porousmediumflow/fv1pconstantspatialparams.hh"
#include <dumux/porousmediumflow/fv1pconstantspatialparams.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
using FVSpatialParamsOnePConstant
[[deprecated("Use FVSpatialParamsOnePConstant in dumux/porousmediumflow/fv1pconstantspatialparams.hh instead!")]]
= FVPorousMediumOnePConstantSpatialParams<GridGeometry, Scalar>;

} // end namespace Dumux

#endif
