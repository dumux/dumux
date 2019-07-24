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
 * \brief The base class for spatial parameters of linear elastic geomechanical problems
 */
#ifndef DUMUX_GEOMECHANICS_ELASTIC_FV_SPATIAL_PARAMS_HH
#define DUMUX_GEOMECHANICS_ELASTIC_FV_SPATIAL_PARAMS_HH

#warning "This header is deprecated and will be removed after release 3.1. Use spatialparams/elastic.hh"

#include "elastic.hh"

namespace Dumux {

// Deprecated class signature
template<class Scalar, class GridGeometry, class Implementation>
using FVSpatialParamsElastic [[deprecated("Use SpatialParamsElastic instead. Will be removed after 3.1!")]]
    = SpatialParamsElastic<Scalar, GridGeometry, Implementation>;

} // end namespace Dumux

#endif
