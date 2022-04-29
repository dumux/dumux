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
#ifndef DUMUX_PYTHON_COMMON_FVSPATIALPARAMS_HH
#define DUMUX_PYTHON_COMMON_FVSPATIALPARAMS_HH

#include <memory>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dune/python/pybind11/pybind11.h>

namespace Dumux::Python {

/*!
 * \ingroup Common
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters used with finite-volume schemes.
 */
template<class GridGeometry_>
class FVSpatialParams
{
public:
    using GridGeometry = GridGeometry_;
    using GridView = typename GridGeometry::GridView;
    using Scalar = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethods::box;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;

    FVSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                    pybind11::object pySpatialParameters)
    : gridGeometry_(gridGeometry)
    , pySpatialParameters_(pySpatialParameters)
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
        if (pybind11::hasattr(pySpatialParameters_, "extrusionFactor"))
            return pySpatialParameters_.attr("extrusionFactor")(element, scv, elemSol).template cast<Scalar>();
        else if (pybind11::hasattr(pySpatialParameters_, "extrusionFactorAtPos"))
            return pySpatialParameters_.attr("extrusionFactorAtPos")(scv.dofPosition()).template cast<Scalar>();
        else
            return 1.0;
    }

    /*!
     * \brief Return the temperature in the given sub-control volume.
     */
    template<class ElementSolution>
    Scalar temperature(const Element& element,
                       const SubControlVolume& scv,
                       const ElementSolution& elemSol) const
    {
        if (pybind11::hasattr(pySpatialParameters_, "temperature"))
            return pySpatialParameters_.attr("temperature")(element, scv, elemSol).template cast<Scalar>();
        else if (pybind11::hasattr(pySpatialParameters_, "temperatureAtPos"))
            return pySpatialParameters_.attr("temperatureAtPos")(scv.dofPosition()).template cast<Scalar>();
        else
            return 283.15;
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
    pybind11::object pySpatialParameters_;
    GravityVector gravity_;
};

// Python wrapper for the above FVProblem C++ class
template<class SpatialParams, class... options>
void registerFVSpatialParams(pybind11::handle scope, pybind11::class_<SpatialParams, options...> cls)
{
    using pybind11::operator""_a;

    using GridGeometry = typename SpatialParams::GridGeometry;

    cls.def(pybind11::init([](std::shared_ptr<const GridGeometry> gridGeometry, pybind11::object p){
        return std::make_shared<SpatialParams>(gridGeometry, p);
    }));

    cls.def("extrusionFactor", &SpatialParams::template extrusionFactor<decltype(std::ignore)>);
    cls.def("temperature", &SpatialParams::template temperature<decltype(std::ignore)>);
    cls.def("gravity", &SpatialParams::gravity);
    cls.def_property_readonly("gridGeometry", &SpatialParams::gridGeometry);
}

} // namespace Dumux::Python

#endif
