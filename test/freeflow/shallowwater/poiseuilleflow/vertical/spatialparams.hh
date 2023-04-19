// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters for the Poiseuille flow problem.
 */
#ifndef DUMUX_POISEUILLE_FLOW_SPATIAL_PARAMETERS_HH
#define DUMUX_POISEUILLE_FLOW_SPATIAL_PARAMETERS_HH

#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/freeflow/spatialparams.hh>

#include <dumux/material/fluidmatrixinteractions/frictionlaws/frictionlaw.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/viscousnoslip.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters class for the Poiseuille flow test.
 *
 */
template<class GridGeometry, class Scalar, class VolumeVariables>
class PoiseuilleFlowSpatialParams
: public FreeFlowSpatialParams<GridGeometry, Scalar,
                         PoiseuilleFlowSpatialParams<GridGeometry, Scalar, VolumeVariables>>
{
    using ThisType = PoiseuilleFlowSpatialParams<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    PoiseuilleFlowSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        gravity_ = getParam<Scalar>("Problem.Gravity");
        bedSlope_ = getParam<Scalar>("Problem.BedSlope");
        channelLength_ = this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0];
        frictionLaw_ = std::make_unique<FrictionLawViscousNoSlip<VolumeVariables>>();
    }

    Scalar gravity(const GlobalPosition& globalPos) const
    { return gravity_; }

    Scalar bedSurface(const Element& element, const SubControlVolume& scv) const
    { return bedSlope_*(channelLength_ - scv.center()[0]); }

    const FrictionLaw<VolumeVariables>& frictionLaw(const Element& element, const SubControlVolume& scv) const
    { return *frictionLaw_; }

private:
    Scalar gravity_;
    Scalar bedSlope_;
    Scalar channelLength_;

    std::unique_ptr<FrictionLaw<VolumeVariables>> frictionLaw_;
};

} // end namespace Dumux

#endif
