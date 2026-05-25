// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTIC_INCOMPRESSIBLE_BLOCK_SPATIAL_PARAMS_HH
#define DUMUX_HYPERELASTIC_INCOMPRESSIBLE_BLOCK_SPATIAL_PARAMS_HH

#include <dumux/common/fvspatialparams.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

// Spatial params that read the Lamé parameters λ and µ directly from the input.
// This avoids E/ν round-tripping for nearly-incompressible benchmark materials.
template<class GridGeometry, class Scalar>
class IncompressibleBlockSpatialParams
: public FVSpatialParams<GridGeometry, Scalar, IncompressibleBlockSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = FVSpatialParams<GridGeometry, Scalar,
                                       IncompressibleBlockSpatialParams<GridGeometry, Scalar>>;
public:
    IncompressibleBlockSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , lambda_(getParam<Scalar>("SpatialParams.Lambda"))
    , mu_(getParam<Scalar>("SpatialParams.Mu"))
    {}

    Scalar shearModulus() const { return mu_; }
    Scalar bulkModulus() const { return lambda_ + 2.0*mu_/3.0; }
    Scalar firstLameParameter() const { return lambda_; }


private:
    Scalar lambda_, mu_;
};

} // end namespace Dumux
#endif
