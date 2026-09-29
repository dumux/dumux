// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters of Wooding's V-catchment for the shallow water model
 */
#ifndef DUMUX_TEST_FREEFLOW_SHALLOWWATER_LONGWAVE_WOODING_SPATIALPARAMS_LOCALINERTIAL_HH
#define DUMUX_TEST_FREEFLOW_SHALLOWWATER_LONGWAVE_WOODING_SPATIALPARAMS_LOCALINERTIAL_HH

#include <memory>

#include <dumux/common/parameters.hh>
#include <dumux/freeflow/spatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/manning.hh>

#include "catchment.hh"

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters of Wooding's V-catchment for the shallow water model
 *
 * The roughness is piecewise constant, so one friction law per material suffices.
 */
template<class GridGeometry, class Scalar, class VolumeVariables>
class WoodingLocalInertialSpatialParams
: public FreeFlowSpatialParams<GridGeometry, Scalar,
                               WoodingLocalInertialSpatialParams<GridGeometry, Scalar, VolumeVariables>>
{
    using ThisType = WoodingLocalInertialSpatialParams<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, ThisType>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    WoodingLocalInertialSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , gravity_(getParam<Scalar>("Problem.Gravity", 9.81))
    , plane_(gravity_, getParam<Scalar>("Problem.PlaneManningN"))
    , channel_(gravity_, getParam<Scalar>("Problem.ChannelManningN"))
    {}

    const FrictionLaw<VolumeVariables>& frictionLaw(const Element& element,
                                                    const SubControlVolume& scv) const
    { return catchment_.inChannel(scv.center()) ? channel_ : plane_; }

    Scalar gravity(const GlobalPosition& globalPos) const
    { return gravity_; }

    Scalar bedSurface(const Element& element, const SubControlVolume& scv) const
    { return catchment_.bedElevation(scv.center()); }

    const Wooding::Catchment<Scalar>& catchment() const
    { return catchment_; }

private:
    Scalar gravity_;
    FrictionLawManning<VolumeVariables> plane_, channel_;
    Wooding::Catchment<Scalar> catchment_;
};

} // end namespace Dumux

#endif
