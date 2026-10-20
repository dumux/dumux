// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters for the manufactured solution test
 */
#ifndef DUMUX_SHALLOWWATER_MANUFACTURED_SPATIAL_PARAMETERS_HH
#define DUMUX_SHALLOWWATER_MANUFACTURED_SPATIAL_PARAMETERS_HH

#include <memory>

#include <dumux/common/parameters.hh>
#include <dumux/freeflow/spatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/frictionlaw.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/nofriction.hh>

#include "exactsolution.hh"

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters for the manufactured solution test: the bed of the manufactured solution, no friction
 */
template<class GridGeometry, class Scalar, class VolumeVariables>
class ManufacturedSolutionSpatialParams
: public FreeFlowSpatialParams<GridGeometry, Scalar, ManufacturedSolutionSpatialParams<GridGeometry, Scalar, VolumeVariables>>
{
    using ThisType = ManufacturedSolutionSpatialParams<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    ManufacturedSolutionSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , gravity_(getParam<Scalar>("Problem.Gravity"))
    , frictionLaw_(std::make_unique<FrictionLawNoFriction<VolumeVariables>>())
    {}

    Scalar gravity(const GlobalPosition& globalPos) const
    { return gravity_; }

    const FrictionLaw<VolumeVariables>& frictionLaw(const Element& element,
                                                    const SubControlVolume& scv) const
    { return *frictionLaw_; }

    Scalar bedSurface(const Element& element,
                      const SubControlVolume& scv) const
    { return bedSurfaceAtPos(scv.center()); }

    Scalar bedSurfaceAtPos(const GlobalPosition& globalPos) const
    { return ShallowWaterTest::ManufacturedSolution<Scalar>::state(globalPos)[3]; }

private:
    Scalar gravity_;
    std::unique_ptr<FrictionLaw<VolumeVariables>> frictionLaw_;
};

} // end namespace Dumux

#endif
