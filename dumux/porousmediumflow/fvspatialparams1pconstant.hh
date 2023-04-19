// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PorousmediumflowModels
 * \ingroup SpatialParameters
 * \brief A spatial params implementation for 1p problem with constant properties
 */
#ifndef DUMUX_POROUS_MEDIUM_FLOW_FV_SPATIAL_PARAMS_ONEP_CONSTANT_HH
#define DUMUX_POROUS_MEDIUM_FLOW_FV_SPATIAL_PARAMS_ONEP_CONSTANT_HH

#include <dumux/common/parameters.hh>
#include "fvspatialparams1p.hh"

namespace Dumux {

/*!
 * \ingroup PorousmediumflowModels
 * \ingroup SpatialParameters
 * \brief A spatial params implementation for 1p problem with constant properties
 */
template<class GridGeometry, class Scalar>
class FVPorousMediumFlowSpatialParamsOnePConstant
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, FVPorousMediumFlowSpatialParamsOnePConstant<GridGeometry, Scalar>>
{
    using ThisType = FVPorousMediumFlowSpatialParamsOnePConstant<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;

    FVPorousMediumFlowSpatialParamsOnePConstant(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , porosity_(getParam<Scalar>("SpatialParams.Porosity"))
    , permeability_(getParam<Scalar>("SpatialParams.Permeability"))
    , temperature_(getParam<Scalar>(
        "SpatialParams.Temperature",
        ParentType::temperatureAtPos(GlobalPosition(0.0))
    ))
    {}

    /*!
     * \brief The (intrinsic) permeability \f$[m^2]\f$
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    /*!
     * \brief The porosity \f$[-]\f$
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    /*!
     * \brief The temperature \f$[K]\f$
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return temperature_; }

private:
    const Scalar porosity_;
    const Scalar permeability_;
    const Scalar temperature_;
};

} // end namespace Dumux

#endif
