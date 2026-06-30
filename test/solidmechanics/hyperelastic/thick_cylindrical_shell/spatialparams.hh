// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTIC_THICK_CYLINDRICAL_SHELL_SPATIAL_PARAMS_HH
#define DUMUX_HYPERELASTIC_THICK_CYLINDRICAL_SHELL_SPATIAL_PARAMS_HH

#include <dumux/common/fvspatialparams.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

// Spatial params for the thick cylindrical shell benchmark.
// The benchmark (Reese et al. 2000 / Elguedj et al. 2008) specifies the shear
// modulus µ and the Lamé parameter Λ directly:
//   µ = 6000 N/mm², Λ = 240000 N/mm²  (=> ν ≈ 0.488, E ≈ 17.85 GPa).
// We therefore read µ and Λ straight from the input rather than E and ν.
template<class GridGeometry, class Scalar>
class ThickCylindricalShellSpatialParams
: public FVSpatialParams<GridGeometry, Scalar,
                         ThickCylindricalShellSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = FVSpatialParams<GridGeometry, Scalar,
                                       ThickCylindricalShellSpatialParams<GridGeometry, Scalar>>;
public:
    ThickCylindricalShellSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , mu_(getParam<Scalar>("SpatialParams.ShearModulus"))
    , lambda_(getParam<Scalar>("SpatialParams.FirstLameParameter"))
    {}

    Scalar shearModulus() const { return mu_; }
    Scalar firstLameParameter() const { return lambda_; }
    Scalar bulkModulus() const { return lambda_ + 2.0*mu_/3.0; }

private:
    Scalar mu_, lambda_;
};

} // end namespace Dumux
#endif
