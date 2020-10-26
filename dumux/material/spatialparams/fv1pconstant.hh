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
 * \ingroup SpatialParameters
 * \brief A spatial params implementation for 1p problem with constant properties
 */
#ifndef DUMUX_FV_CONSTANT_SPATIAL_PARAMS_ONE_P_HH
#define DUMUX_FV_CONSTANT_SPATIAL_PARAMS_ONE_P_HH

#include <dumux/common/parameters.hh>
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

/*!
 * \ingroup SpatialParameters
 * \brief A spatial params implementation for 1p problem with constant properties
 */
template<class GridGeometry, class Scalar>
class FVSpatialParamsOnePConstant
: public FVSpatialParamsOneP<GridGeometry, Scalar, FVSpatialParamsOnePConstant<GridGeometry, Scalar>>
{
    using ThisType = FVSpatialParamsOnePConstant<GridGeometry, Scalar>;
    using ParentType = FVSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;

    FVSpatialParamsOnePConstant(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , porosity_(getParam<Scalar>("SpatialParams.Porosity"))
    , permeability_(getParam<Scalar>("SpatialParams.Permeability"))
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

private:
    const Scalar porosity_;
    const Scalar permeability_;
};

} // end namespace Dumux

#endif
