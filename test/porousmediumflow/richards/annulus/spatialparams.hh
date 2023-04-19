// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_RICHARDS_ANNULUS_SPATIAL_PARAMETERS_HH
#define DUMUX_RICHARDS_ANNULUS_SPATIAL_PARAMETERS_HH

#include <iostream>
#include <dumux/common/parameters.hh>

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

namespace Dumux {

/*!
 * \ingroup RichardsTests
 * \brief Spatial parameters for the Richards benchmarks
 */
template<class GridGeometry, class Scalar>
class RichardsAnnulusSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, RichardsAnnulusSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = RichardsAnnulusSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;
    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    RichardsAnnulusSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , paramGroup_("SpatialParams." + getParam<std::string>("SpatialParams.SoilType"))
    , pcKrSwCurve_(paramGroup_)
    , permeability_(getParam<Scalar>(paramGroup_ + ".Permeability"))
    , porosity_(getParam<Scalar>(paramGroup_ + ".Porosity"))
    {
        std::cout << "Using " << getParam<std::string>(paramGroup_ + ".Name") << " parameters" << std::endl;
    }

    /*!
     * \brief Returns the intrinsic permeability tensor [m^2] at a given location
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    /*!
     * \brief Returns the porosity [-] at a given location
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    /*!
     * \brief Returns the fluid-matrix interaction law for the sub-control volume
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    { return makeFluidMatrixInteraction(pcKrSwCurve_); }

    /*!
     * \brief material law for analytic solution
     */
    const PcKrSwCurve& pcKrSwCurve() const { return pcKrSwCurve_; }

    /*!
     * \brief the parameter group selecting the soil type
     */
    const std::string& paramGroup() const { return paramGroup_; }

private:
    std::string paramGroup_;
    const PcKrSwCurve pcKrSwCurve_;
    const Scalar permeability_, porosity_;
};

} // end namespace Dumux

#endif
