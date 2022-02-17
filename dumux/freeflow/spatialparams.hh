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
 * \ingroup FreeflowModels
 * \brief Definition of the spatial parameters for the freeflow problems.
 */
#ifndef DUMUX_FREEFLOW_SPATIAL_PARAMS_HH
#define DUMUX_FREEFLOW_SPATIAL_PARAMS_HH

#include <dumux/common/fvspatialparams.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Definition of the spatial parameters for the freeflow problems.
 */
template<class GridGeometry, class Scalar, class Implementation>
class FreeFlowSpatialParams
: public FVSpatialParams<GridGeometry, Scalar, Implementation>
{
    using ParentType = FVSpatialParams<GridGeometry, Scalar, Implementation>;

public:
    FreeFlowSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}
};

/*!
 * \ingroup FreeflowModels
 * \brief Definition of the spatial parameters for the freeflow problems.
 */
template<class GridGeometry, class Scalar>
class FreeFlowDefaultSpatialParams
: public FreeFlowSpatialParams<GridGeometry, Scalar, FreeFlowDefaultSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, FreeFlowDefaultSpatialParams<GridGeometry, Scalar>>;
public:
    using ParentType::ParentType;
};

} // end namespace Dumux

#endif
