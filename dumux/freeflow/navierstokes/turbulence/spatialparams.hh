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
 * \ingroup RANSSpatialParameters
 * \brief The base class for spatial parameters of the RANS models
 * using a fully implicit discretization method.
 */

#ifndef DUMUX_FREEFLOW_TURBULENCE_SPATIAL_PARAMS_HH
#define DUMUX_FREEFLOW_TURBULENCE_SPATIAL_PARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

namespace Dumux {
/*!
 * \ingroup RANSSpatialParameters
 * \brief Definition of the spatial parameters of the RANS models
 */
template<class GridGeometry>
class RANSSpatialParams
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Scalar = typename GlobalPosition::value_type;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int dim = GlobalPosition::dimension;

public:
    RANSSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    {
        fixedFlowDirectionAxis_ = getParam<int>("RANS.FlowDirectionAxis", 0);
        fixedWallNormalAxis_ = getParam<int>("RANS.WallNormalAxis", 1);
    }

    void setWallDistance(const std::vector<Scalar>& wallDistance)
    { wallDistance_ = wallDistance; }

    Scalar wallDistance(const int elementIdx) const
    { return wallDistance_[elementIdx]; }

    template <class WallData>
    void setWallData(const WallData& wallData, const GridGeometry &gridGeometry)
    {
        wallNormalAxis_.resize(wallData.size());
        wallElementIdx_.resize(wallData.size());
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            unsigned int elementIdx = gridGeometry.elementMapper().index(element);
            wallElementIdx_[elementIdx] = wallData[elementIdx].eIdx;
            if ( ! (hasParam("RANS.WallNormalAxis")) )
            {
                GlobalPosition wallOuterNormal = wallData[elementIdx].scvfOuterNormal;
                if constexpr (dim == 2) // 2D
                    wallNormalAxis_[elementIdx] = (wallOuterNormal[0] == 1) ? 0 : 1;
                else // 3D
                    wallNormalAxis_[elementIdx] = (wallOuterNormal[0] == 1) ? 0 : ((wallOuterNormal[1] == 1) ? 1 : 2);
            }
            else
                wallNormalAxis_[elementIdx] = fixedWallNormalAxis_;
        }
    }

private:
    std::vector<Scalar> wallDistance_;
    std::vector<unsigned int> wallNormalAxis_;
    std::vector<unsigned int> flowDirectionAxis_;
    std::vector<unsigned int> wallElementIdx_;

    int fixedFlowDirectionAxis_;
    int fixedWallNormalAxis_;

};

}

#endif
