// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Spatial parameters for the dimensionless Boussinesq dissolution test.
 *
 * Permeability and porosity are both 1 (dimensionless).
 * Gravity is set to unit magnitude pointing downward, overriding the base
 * class default of 9.81 m/s².
 */
#ifndef DUMUX_BOUSSINESQ_SPATIAL_PARAMS_HH
#define DUMUX_BOUSSINESQ_SPATIAL_PARAMS_HH

#include <dune/common/fvector.hh>

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class BoussinesqSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                             BoussinesqSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                           BoussinesqSpatialParams<GridGeometry, Scalar>>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using GravityVector  = Dune::FieldVector<Scalar, dimWorld>;

public:
    using PermeabilityType = Scalar;

    BoussinesqSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        gravity_ = Scalar(0);
        gravity_[dimWorld-1] = -1.0; // dimensionless unit gravity pointing down
    }

    PermeabilityType permeabilityAtPos(const GlobalPosition&) const
    { return 1.0; }

    Scalar porosityAtPos(const GlobalPosition&) const
    { return 1.0; }

    //! Override base-class gravity with the dimensionless unit vector
    const GravityVector& gravity(const GlobalPosition&) const
    { return gravity_; }

private:
    GravityVector gravity_;
};

} // namespace Dumux

#endif