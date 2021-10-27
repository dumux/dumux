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
 * \ingroup PorousMediumFlow
 * \ingroup SpatialParameters
 * \brief A spatial params implementation for 1p problem with constant properties
 */
#ifndef DUMUX_POROUS_MEDIUM_FLOW_FV_SPATIAL_PARAMS_ONEP_CONSTANT_HH
#define DUMUX_POROUS_MEDIUM_FLOW_FV_SPATIAL_PARAMS_ONEP_CONSTANT_HH

#include <dumux/common/parameters.hh>
#include "fvspatialparams1p.hh"

namespace Dumux {

/*!
 * \ingroup SpatialParameters
 * \brief A spatial params implementation for 1p problem with constant properties
 */
template<class GridGeometry, class Scalar>
class FVPorousMediumSpatialParamsOnePConstant
: public FVPorousMediumSpatialParamsOneP<GridGeometry, Scalar, FVPorousMediumSpatialParamsOnePConstant<GridGeometry, Scalar>>
{
    using ThisType = FVPorousMediumSpatialParamsOnePConstant<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;

    FVPorousMediumSpatialParamsOnePConstant(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , porosity_(getParam<Scalar>("SpatialParams.Porosity"))
    , permeability_(getParam<Scalar>("SpatialParams.Permeability"))
    , temperature_(getParam<Scalar>("SpatialParams.Temperature"))
    {}

    /*!
     * \brief The (intrinsic) permeability \f$[m^2]\f$
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    /*!
     * \brief The porosity \f$[-]\f$
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    /*!
     * \brief The temperature \f$[K]\f$
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return temperature_; }

private:
    const Scalar porosity_;
    const Scalar permeability_;
    const Scalar temperature_;
};

} // end namespace Dumux

#endif
