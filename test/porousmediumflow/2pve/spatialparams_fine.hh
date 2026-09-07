// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup VETest
 * \brief The fine-level spatial params for the vertical equilibrium Darcy test.
 */

#ifndef DUMUX_TEST_TWOPVE_SPATIALPARAMS_FINE_HH
#define DUMUX_TEST_TWOPVE_SPATIALPARAMS_FINE_HH

#include <memory>

#include <dune/common/fvector.hh>

#include <dumux/common/parameters.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class TwoPTestFineSpatialParams
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;

public:
    using PermeabilityType = Scalar;

    TwoPTestFineSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
        : gridGeometry_(gridGeometry),
          permeability_ (getParam<Scalar>("SpatialParams.Permeability")),
          porosity_ (getParam<Scalar>("SpatialParams.Porosity")),
          gravity_(0.0)
    {
        if (getParam<bool>("Problem.EnableGravity"))
            gravity_[dimWorld-1] = -9.81;
    }


    /*!
     * \brief Returns the intrinsic permeability tensor \f$[m^2]\f$ for a specific fine-level element
     *
     * \param element fine-level element
     */
    decltype(auto) permeabilityAtElement(const Element& element) const
    {
        if constexpr(dimWorld==2)
        {
            const auto pos = element.geometry().center();
            const auto relPos = pos - gridGeometry_->bBoxMin();

            bool isInUpperLens = relPos[1] > 45.0 && relPos[0] > 70.0 && relPos[0] < 85.0;
            if(isInUpperLens)
                return permeability_*0.00001;
            else
                return permeability_;
        }
        else if constexpr(dimWorld==3)
        {
            const auto pos = element.geometry().center();
            const auto relPos = pos - gridGeometry_->bBoxMin();

            bool isInLens = relPos[0]>20.0 && relPos[0]<80 && relPos[1]>35.0 && relPos[1]<65.0 && relPos[dimWorld-1]>15.0;
            if(isInLens)
                return permeability_*0.00001;
            else
                return permeability_;
        }
    }


    /*!
     * \brief Returns the porosity for a specific fine-level element
     *
     * \param element fine-level element
     */
    Scalar porosityAtElement(const Element& element) const
    {
        return porosity_;
    }


    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * The default behaviour is a constant gravity vector;
     * if the <tt>Problem.EnableGravity</tt> parameter is true,
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$,
     * else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$.
     *
     * \param pos the spatial position at which to evaluate the gravity vector
     */
    const GravityVector& gravity(const GlobalPosition& pos) const
    { return gravity_; }

    //! The finite volume grid geometry
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
    PermeabilityType permeability_;
    Scalar porosity_;
    GravityVector gravity_;
};

} // end namespace Dumux

#endif
