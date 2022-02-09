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
 * \ingroup PoreNetworkModels
 * \ingroup SpatialParameters
 * \brief The spatial parameters for single-phase pore-network models.
 */
#ifndef DUMUX_PNM_1P_SPATIAL_PARAMS_HH
#define DUMUX_PNM_1P_SPATIAL_PARAMS_HH

#include <dumux/porenetwork/common/spatialparams.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PoreNetworkModels
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
 * \ingroup PoreNetworkModels
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
