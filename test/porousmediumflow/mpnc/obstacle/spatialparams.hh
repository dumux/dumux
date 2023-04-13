// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MPNCTests
 * \brief The spatial parameters for the ObstacleProblem.
 */

#ifndef DUMUX_OBSTACLE_SPATIAL_PARAMS_HH
#define DUMUX_OBSTACLE_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>
#include <dumux/material/fluidmatrixinteractions/2p/smoothedlinearlaw.hh>
#include <dumux/material/fluidmatrixinteractions/mp/mpadapter.hh>

namespace Dumux {

/**
 * \ingroup MPNCTests
 * \brief Definition of the spatial params properties for the obstacle problem.
 *
 */
template<class GridGeometry, class Scalar>
class ObstacleSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                       ObstacleSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                                     ObstacleSpatialParams<GridGeometry, Scalar>>;

    enum {dimWorld=GridView::dimensionworld};
    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    using PcKrSwCurve = FluidMatrix::SmoothedLinearLaw<Scalar>;

public:
    //! Export the type used for the permeability
    using PermeabilityType = Scalar;


    ObstacleSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwCurve_("SpatialParams") // initialize the material law
    , coarseK_(1e-12) // intrinsic permeability
    , fineK_(1e-15) // intrinsic permeability
    , porosity_(0.3)
    , temperature_(getParam<Scalar>("SpatialParams.Temperature"))
    {}

    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        if (isFineMaterial_(scv.dofPosition()))
            return fineK_;
        else
            return coarseK_;
    }

    /*!
     * \brief Defines the porosity \f$[-]\f$ of the soil
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    /*!
     * \brief Defines the temperature \f$[K]\f$ at the given position
     * \param globalPos The global position
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return temperature_; }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition &globalPos) const
    {
        return makeFluidMatrixInteraction(FluidMatrix::MPAdapter(pcKrSwCurve_));
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The global position
     * \return The wetting phase index
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        return FluidSystem::phase0Idx;
    }

private:
    /*!
     * \brief Returns whether a given global position is in the
     *        fine-permeability region or not.
     */
    static bool isFineMaterial_(const GlobalPosition &pos)
    {
        return
            10 - eps_ <= pos[0] && pos[0] <= 20 + eps_ &&
            0 - eps_ <= pos[1] && pos[1] <= 35 + eps_;
    }

    const PcKrSwCurve pcKrSwCurve_;
    const Scalar coarseK_;
    const Scalar fineK_;
    const Scalar porosity_;
    const Scalar temperature_;
    static constexpr Scalar eps_ = 1e-6;
};

} // end namespace Dumux

#endif
