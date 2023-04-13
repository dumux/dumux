// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
