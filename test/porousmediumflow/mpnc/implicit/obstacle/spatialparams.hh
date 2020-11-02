// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup MPNCTests
 * \brief The spatial parameters for the ObstacleProblem.
 */

#ifndef DUMUX_OBSTACLE_SPATIAL_PARAMS_HH
#define DUMUX_OBSTACLE_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv.hh>
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
: public FVSpatialParams<GridGeometry, Scalar,
                         ObstacleSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParams<GridGeometry, Scalar,
                                       ObstacleSpatialParams<GridGeometry, Scalar>>;

    enum {dimWorld=GridView::dimensionworld};
    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    using PcKrSwCurve = FluidMatrix::SmoothedLinearLaw<Scalar>;
    using MPAdapter = Dumux::FluidMatrix::MPAdapter<PcKrSwCurve, 2>;

public:
    //! Export the type used for the permeability
    using PermeabilityType = Scalar;


    ObstacleSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwCurve_("SpatialParams") // initialize the material law
    , coarseK_(1e-12) // intrinsic permeability
    , fineK_(1e-15) // intrinsic permeability
    , porosity_(0.3)
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
     * \brief Returns the fluid-matrix interaction law at a given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition &globalPos) const
    {
        return makeFluidMatrixInteraction(MPAdapter(pcKrSwCurve_));
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
    static constexpr Scalar eps_ = 1e-6;
};

} // end namespace Dumux

#endif
