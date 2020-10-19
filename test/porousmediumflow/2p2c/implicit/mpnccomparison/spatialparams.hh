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
 * \ingroup TwoPTwoCTests
 * \brief The spatial parameters for the 2p2c mpnc comparison problem.
 */

#ifndef DUMUX_MPNC_COMPARISON_SPATIAL_PARAMS_HH
#define DUMUX_MPNC_COMPARISON_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

namespace Dumux {

/**
 * \ingroup TwoPTwoCTests
 * \brief The spatial parameters for the 2p2c mpnc comparison problem.
 */
template<class GridGeometry, class Scalar>
class TwoPTwoCComparisonSpatialParams
: public FVSpatialParams<GridGeometry, Scalar,
                         TwoPTwoCComparisonSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParams<GridGeometry, Scalar,
                                       TwoPTwoCComparisonSpatialParams<GridGeometry, Scalar>>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    enum {dimWorld=GridView::dimensionworld};

    using PcKrSwCurve =FluidMatrix::BrooksCoreyDefault<Scalar>;

public:
    //! Export permeability type
    using PermeabilityType = Scalar;

    TwoPTwoCComparisonSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , finePcKrSwCurve_("SpatialParams.FineMaterial")
    , coarsePcKrSwCurve_("SpatialParams.CoarseMaterial")
    {
        // intrinsic permeabilities
        coarseK_ = 1e-12;
        fineK_ = 1e-15;

        // the porosity
        porosity_ = 0.3;
    }

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
     * \param globalPos The global position of the sub-control volume.
     * \return the material parameters object
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_;
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param globalPos The global position of the sub-control volume.
     * \return The material parameters object
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return makeFluidMatrixInteraction(finePcKrSwCurve_);
        return makeFluidMatrixInteraction(coarsePcKrSwCurve_);
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The global position of the sub-control volume.
     * \return The wetting phase index
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::H2OIdx; }

private:
    /*!
     * \brief Returns whether a given global position is in the
     *        fine-permeability region or not.
     * \param pos The global position of the sub-control volume.
     */
    static bool isFineMaterial_(const GlobalPosition &pos)
    {
        return
            30 - eps_ <= pos[0] && pos[0] <= 50 + eps_ &&
            20 - eps_ <= pos[1] && pos[1] <= 40 + eps_;
    }

    Scalar coarseK_;
    Scalar fineK_;
    Scalar porosity_;

    const PcKrSwCurve finePcKrSwCurve_;
    const PcKrSwCurve coarsePcKrSwCurve_;
    static constexpr Scalar eps_ = 1e-6;
};

} // end namespace Dumux

#endif
