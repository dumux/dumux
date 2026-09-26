// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVETests
 * \brief The spatial parameters for the radial injection into a homogeneous confined aquifer.
 */

#ifndef DUMUX_TEST_TWOPVE_RADIAL_INJECTION_SPATIAL_PARAMS_HH
#define DUMUX_TEST_TWOPVE_RADIAL_INJECTION_SPATIAL_PARAMS_HH

#include <memory>

#include <dumux/common/parameters.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>
#include <dumux/porousmediumflow/2pve/columnmapping.hh>
#include <dumux/porousmediumflow/2pve/spatialparams.hh>

namespace Dumux {

/*!
 * \ingroup TwoPVETests
 * \brief The fine-level spatial parameters of a homogeneous aquifer
 */
template<class GridGeometry, class Scalar>
class TwoPVERadialInjectionFineSpatialParams
{
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

public:
    TwoPVERadialInjectionFineSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_(gridGeometry)
    , permeability_(getParam<Scalar>("SpatialParams.Permeability"))
    , porosity_(getParam<Scalar>("SpatialParams.Porosity"))
    {}

    Scalar permeabilityAtElement(const Element& element) const
    { return permeability_; }

    Scalar porosityAtElement(const Element& element) const
    { return porosity_; }

    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
    Scalar permeability_;
    Scalar porosity_;
};

/*!
 * \ingroup TwoPVETests
 * \brief The coarse-level spatial parameters of a homogeneous aquifer
 *
 * A small Brooks-Corey entry pressure yields a thin capillary fringe,
 * which approximates the sharp interface of the similarity solution.
 */
template<class GridGeometry, class Scalar>
class TwoPVERadialInjectionSpatialParams
: public TwoPVESpatialParams<GridGeometry, Scalar,
                             TwoPVERadialInjectionFineSpatialParams<GridGeometry, Scalar>,
                             TwoPVERadialInjectionSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = TwoPVERadialInjectionSpatialParams<GridGeometry, Scalar>;
    using ParentType = TwoPVESpatialParams<GridGeometry, Scalar, TwoPVERadialInjectionFineSpatialParams<GridGeometry, Scalar>, ThisType>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using PcKrSwCurve = FluidMatrix::BrooksCoreyNoReg<Scalar>;

public:
    using SpatialParamsFine = TwoPVERadialInjectionFineSpatialParams<GridGeometry, Scalar>;

    TwoPVERadialInjectionSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                                       const TwoPVEColumnMapping<GridGeometry, Scalar>& columnMapping,
                                       std::shared_ptr<const SpatialParamsFine> spatialParamsFine,
                                       Scalar fineCellHeight)
    : ParentType(gridGeometry, columnMapping, spatialParamsFine, fineCellHeight)
    , pcKrSwCurve_("SpatialParams")
    {}

    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    { return makeFluidMatrixInteraction(pcKrSwCurve_); }

    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 358.15; }

private:
    const PcKrSwCurve pcKrSwCurve_;
};

} // end namespace Dumux

#endif
