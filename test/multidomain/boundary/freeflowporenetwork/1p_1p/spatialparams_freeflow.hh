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
 * \ingroup BoundaryTests
 * \brief Definition of the spatial parameters for the pore network coupled channel problem.
 */

#ifndef DUMUX_CHANNEL_SPATIAL_PARAMS_HH
#define DUMUX_CHANNEL_SPATIAL_PARAMS_HH

#include <dumux/freeflow/spatialparams.hh>

namespace Dumux {

/*!
 * \ingroup BoundaryTests
 * \brief Definition of the spatial parameters for the pseudo 3D channel freeflow problems.
 */
template<class GridGeometry, class Scalar>
class ChannelSpatialParams
: public FreeFlowSpatialParams<GridGeometry, Scalar, ChannelSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, ChannelSpatialParams<GridGeometry, Scalar>>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    ChannelSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        singleThroatTest_ = getParam<bool>("Problem.SingleThroatTest", true);
        enablePseudoThreeDWallFriction_ = !singleThroatTest_;
        extrusionFactor_ = enablePseudoThreeDWallFriction_ ? getParam<Scalar>("FreeFlow.Problem.Height") : 1.0;
    }

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    { return extrusionFactor_; }

    /*!
     * \brief Return how much the domain is extruded at a given position.
     */
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    { return extrusionFactor_; }

private:
    bool singleThroatTest_;
    bool enablePseudoThreeDWallFriction_;
    Scalar extrusionFactor_;
};

} // end namespace Dumux

#endif
