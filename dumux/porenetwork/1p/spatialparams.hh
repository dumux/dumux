// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMOnePModel
 * \ingroup SpatialParameters
 * \brief The spatial parameters for single-phase pore-network models.
 */
#ifndef DUMUX_PNM_1P_SPATIAL_PARAMS_HH
#define DUMUX_PNM_1P_SPATIAL_PARAMS_HH

#include <dumux/porenetwork/common/spatialparams.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup SpatialParameters
 * \ingroup PNMOnePModel
 * \brief The base class for spatial parameters for single-phase pore-network models.
 */
template<class GridGeometry, class Scalar, class Implementation>
class OnePSpatialParams
: public SpatialParams<GridGeometry, Scalar, Implementation>
{
    using ParentType = SpatialParams<GridGeometry, Scalar, Implementation>;
public:
    using ParentType::ParentType;
};

/*!
 * \ingroup PNMOnePModel
 * \ingroup SpatialParameters
 * \brief The default class for spatial parameters for single-phase pore-network models.
 * \note We have this layer for consistency with the two-phase pore-network models. Also, we
 *       may use this in the feature to define defaults for newly added parameter interfaces.
 */
template<class GridGeometry, class Scalar>
class OnePDefaultSpatialParams
: public OnePSpatialParams<GridGeometry, Scalar, OnePDefaultSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = OnePSpatialParams<GridGeometry, Scalar, OnePDefaultSpatialParams<GridGeometry, Scalar>>;
public:
    using ParentType::ParentType;
};
} // namespace Dumux::PoreNetwork

#endif
