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
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters for the Poiseuille flow problem.
 */
#ifndef DUMUX_POISEUILLE_FLOW_SPATIAL_PARAMETERS_HH
#define DUMUX_POISEUILLE_FLOW_SPATIAL_PARAMETERS_HH

#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/fvspatialparams.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters class for the Poiseuille flow test.
 *
 */
template<class GridGeometry, class Scalar, class VolumeVariables>
class PoiseuilleFlowSpatialParams
: public FVSpatialParams<GridGeometry, Scalar,
                         PoiseuilleFlowSpatialParams<GridGeometry, Scalar, VolumeVariables>>
{
    using ThisType = PoiseuilleFlowSpatialParams<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    PoiseuilleFlowSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        gravity_ = getParam<Scalar>("Problem.Gravity");
        bedSlope_ = getParam<Scalar>("Problem.BedSlope");
        hBoundary_ = getParam<Scalar>("Problem.HBoundary");
    }

    Scalar gravity(const GlobalPosition& globalPos) const
    { return gravity_; }

    Scalar gravity() const
    { return gravity_; }

    //! Define the bed surface
    Scalar bedSurface(const Element& element, const SubControlVolume& scv) const
    {
        const Scalar length = this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0];
        // The bed level sloped downwards from the left boundary to the right boundary
        // The water level boundary condition hBoundary_ is specified at the right boundary
        const Scalar leftBedLevel = - hBoundary_ + length*bedSlope_;
        // todo depends on index e.g. eIdx = scv.elementIndex();
        return leftBedLevel - element.geometry().center()[0] * bedSlope_;
    }

private:
    Scalar gravity_;
    Scalar bedSlope_;
    Scalar hBoundary_;
};

} // end namespace Dumux

#endif
