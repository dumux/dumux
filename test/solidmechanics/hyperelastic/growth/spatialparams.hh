// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTICITY_GROWTH_SPATIAL_PARAMS_HH
#define DUMUX_HYPERELASTICITY_GROWTH_SPATIAL_PARAMS_HH

#include <cmath>
#include <dumux/common/fvspatialparams.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class HyperelasticGrowthSpatialParams
: public FVSpatialParams<GridGeometry, Scalar, HyperelasticGrowthSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = HyperelasticGrowthSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    HyperelasticGrowthSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , eCore_(getParam<Scalar>("SpatialParams.YoungsModulusCore"))
    , eShell_(getParam<Scalar>("SpatialParams.YoungsModulusShell"))
    , nu_(getParam<Scalar>("SpatialParams.PoissonRatio"))
    , rTransition_(getParam<Scalar>("SpatialParams.TransitionRadius"))
    {
        muCore_ = eCore_/(2*(1 + nu_));
        muShell_ = eShell_/(2*(1 + nu_));
        kCore_ = eCore_/(3*(1 - 2*nu_));
        kShell_ = eShell_/(3*(1 - 2*nu_));
    }

    bool isCoreAtPos(const GlobalPosition& globalPos) const
    {
        const auto r = std::hypot(globalPos[0], globalPos[1]);
        return r < rTransition_;
    }

    Scalar shearModulusAtPos(const GlobalPosition& globalPos) const
    { return isCoreAtPos(globalPos) ? muCore_ : muShell_; }

    Scalar bulkModulusAtPos(const GlobalPosition& globalPos) const
    { return isCoreAtPos(globalPos) ? kCore_ : kShell_; }

    Scalar transitionRadius() const
    { return rTransition_; }

private:
    Scalar eCore_, eShell_, nu_;
    Scalar rTransition_;
    Scalar muCore_, muShell_, kCore_, kShell_;
};

} // end namespace Dumux

#endif
