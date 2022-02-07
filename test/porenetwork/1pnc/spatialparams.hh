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
 * \brief Spatial parameters for an isothermal 1p pore-network model
 */
#ifndef DUMUX_COMPOSITIONAL_PNM_SPATIAL_PARAMS_1P_HH
#define DUMUX_COMPOSITIONAL_PNM_SPATIAL_PARAMS_1P_HH

#include <dumux/porenetwork/common/spatialparams.hh>

namespace Dumux::PoreNetwork {

template<class GridGeometry, class Scalar>
class CompositionalSpatialParams
: public SpatialParams<GridGeometry, Scalar, CompositionalSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = SpatialParams<GridGeometry, Scalar, CompositionalSpatialParams<GridGeometry, Scalar>>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
public:
    using PermeabilityType = Scalar;
    using ParentType::ParentType;

    CompositionalSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    { }

    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 273.15 + 10.0; }
};
} // end namespace Dumux::PoreNetwork

#endif
