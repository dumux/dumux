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
 * \brief Basic spatial parameters to be used with finite-volume schemes.
 */
#ifndef DUMUX_COMMON_FV_SPATIAL_PARAMS_HH
#define DUMUX_COMMON_FV_SPATIAL_PARAMS_HH

#include <memory>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters used with finite-volume schemes.
 */
template<class GridGeometry,
         class Scalar,
         class Implementation>
class FVSpatialParams
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;

public:
    FVSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_(gridGeometry)
    , gravity_(0.0)
    {
        if (getParam<bool>("Problem.EnableGravity"))
            gravity_[dimWorld-1] = -9.81;
    }

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    {
        // forward to generic interface
        return asImp_().extrusionFactorAtPos(scv.center());
    }

    /*!
     * \brief Return how much the domain is extruded at a given position.
     */
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    {
        // As a default, i.e. if the user's problem does not overload
        // any extrusion factor method, return 1.0
        return 1.0;
    }

    /*!
     * \brief Return the temperature in the domain at the given position
     * \param globalPos The position in global coordinates where the temperature should be specified.
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { DUNE_THROW(Dune::NotImplemented, "Implementation does not provide temperatureAtPos()"); }

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

protected:
    //! Returns the implementation of the spatial parameters (static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
    GravityVector gravity_;
};

} // namespace Dumux

#endif
