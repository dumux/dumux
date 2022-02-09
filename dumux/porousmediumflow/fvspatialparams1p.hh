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
 * \brief The base class for spatial parameters in single-phase porous-medium-flow problems.
 */
#ifndef DUMUX_POROUS_MEDIUM_FLOW_FV_SPATIAL_PARAMS_ONEP_HH
#define DUMUX_POROUS_MEDIUM_FLOW_FV_SPATIAL_PARAMS_ONEP_HH

#include "fvspatialparams.hh"

namespace Dumux {

/*!
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of single-phase problems.
 */
template<class GridGeometry, class Scalar, class Implementation>
class FVPorousMediumFlowSpatialParamsOneP
: public FVPorousMediumFlowSpatialParams<GridGeometry, Scalar, Implementation>
{
    using ParentType = FVPorousMediumFlowSpatialParams<GridGeometry, Scalar, Implementation>;

public:
    using ParentType::ParentType;
};

} // namespace Dumux

#endif
