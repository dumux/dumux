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
 * \ingroup Common
 * \ingroup SpatialParameters
 * \brief Basic spatial parameters
 */
#ifndef DUMUX_COMMON_FV_SPATIAL_PARAMS_BASE_HH
#define DUMUX_COMMON_FV_SPATIAL_PARAMS_BASE_HH

#include <memory>

#include <dune/common/fvector.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of multi-phase problems
 * using a fully implicit discretization method.
 */
template<class GridGeometry,
         class Scalar = typename GridGeometry::GridView::ctype>
class FVSpatialParamsBase
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;

public:
    FVSpatialParamsBase(std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_(gridGeometry)
    , gravity_(0.0)
    {
        if (getParam<bool>("Problem.EnableGravity"))
            gravity_[dimWorld-1] = -9.81;
    }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * The default behaviour is a constant gravity vector;
     * if the <tt>Problem.EnableGravity</tt> parameter is true,
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$,
     * else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$.
     *
     * \param pos the spatial position at which to evaulate the gravity vector
     */
    const GravityVector& gravity(const GlobalPosition& pos) const
    { return gravity_; }

    //! The finite volume grid geometry
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
    GravityVector gravity_;
};

} // namespace Dumux

#endif
