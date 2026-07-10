// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Spatial parameters for the Boussinesq vorticity test.
 *
 * Provides medium properties (permeability, porosity) and the gravity vector.
 * Fluid properties (μ, ρ₀, β, D) live in the FluidSystem, not here.
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
        // Dimensionless unit gravity pointing downward (last spatial coordinate)
        gravity_           = Scalar(0);
        gravity_[dimWorld-1] = -1.0;
    }

    // ---- permeability (base-class interface) ----
    PermeabilityType permeabilityAtPos(const GlobalPosition&) const
    { return 1.0; }

    // ---- porosity ----
    Scalar porosityAtPos(const GlobalPosition&) const
    { return 1.0; }

    //! Gravity vector [m/s²] — override base-class default (9.81) with unit vector
    const GravityVector& gravity(const GlobalPosition& /*pos*/) const
    { return gravity_; }

private:
    GravityVector gravity_;
};

} // namespace Dumux

#endif
