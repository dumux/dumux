// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief Definition of the spatial parameters for the evaporation constant component problem.
 */
#ifndef DUMUX_EVAPORATION_CONSTANT_COMPONENT_SPATIAL_PARAMS_HH
#define DUMUX_EVAPORATION_CONSTANT_COMPONENT_SPATIAL_PARAMS_HH

#include <dumux/common/math.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/ploteffectivediffusivitymodel.hh>
#include <dumux/io/plotpckrsw.hh>
#include <dumux/io/plotthermalconductivitymodel.hh>
#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTwoCModel
 * \brief Definition of the spatial parameters for the evaporation constant component problem.
 */
template<class GridGeometry, class Scalar>
class EvaporationConstantComponentSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                         EvaporationConstantComponentSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = EvaporationConstantComponentSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using Element = typename GridView::template Codim<0>::Entity;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

public:
    //! Export the type used for the permeability
    using PermeabilityType = Scalar;


    EvaporationConstantComponentSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwCurve_("SpatialParams")
    {
        // intrinsic permeability
        permeability_ = 1e-11;

        // porosity
        porosity_ = 0.3;
    }


    /*!
     * \brief Applies the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param globalPos The global position
     */
    Scalar permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        return permeability_;
    }

    /*!
     * \brief Defines the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_;
    }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     * \param globalPos The global coordinates for the given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        return makeFluidMatrixInteraction(pcKrSwCurve_);
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The position of the center of the element
     * \return The wetting phase index
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }

private:
    Scalar permeability_;
    Scalar porosity_;

    const PcKrSwCurve pcKrSwCurve_;
};

} // end namespace Dumux

#endif
