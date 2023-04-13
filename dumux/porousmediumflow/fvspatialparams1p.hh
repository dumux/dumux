// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
