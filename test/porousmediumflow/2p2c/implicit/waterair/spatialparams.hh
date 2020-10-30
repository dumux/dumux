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
 * \brief Definition of the spatial parameters for the water-air problem.
 */
#ifndef DUMUX_WATER_AIR_SPATIAL_PARAMS_HH
#define DUMUX_WATER_AIR_SPATIAL_PARAMS_HH

#include <dumux/common/math.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/ploteffectivediffusivitymodel.hh>
#include <dumux/io/plotpckrsw.hh>
#include <dumux/io/plotthermalconductivitymodel.hh>
#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTwoCModel
 * \brief Definition of the spatial parameters for the water-air problem.
 */
template<class GridGeometry, class Scalar>
class WaterAirSpatialParams
: public FVSpatialParams<GridGeometry, Scalar,
                         WaterAirSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = WaterAirSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using Element = typename GridView::template Codim<0>::Entity;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

public:
    //! Export the type used for the permeability
    using PermeabilityType = Scalar;


    WaterAirSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwCurve_("SpatialParams")
    {
        layerBottom_ = 22.0;

        // intrinsic permeabilities
        fineK_ = 1e-13;
        coarseK_ = 1e-12;

        // porosities
        finePorosity_ = 0.3;
        coarsePorosity_ = 0.3;

        plotFluidMatrixInteractions_ = getParam<bool>("Output.PlotFluidMatrixInteractions");
    }

    /*!
     * \brief This is called from the problem and creates a gnuplot output
     *        of e.g the pc-Sw curve
     */
    void plotMaterialLaw()
    {
        GnuplotInterface<Scalar> gnuplot(plotFluidMatrixInteractions_);
        gnuplot.setOpenPlotWindow(plotFluidMatrixInteractions_);

        const auto sw = linspace(0.2, 1.0, 1000);

        const auto pc = samplePcSw(pcKrSwCurve_, sw);
        Gnuplot::addPcSw(gnuplot, sw, pc, "pc-Sw", "w lp");
        gnuplot.setOption("set xrange [0:1]");
        gnuplot.setOption("set label \"residual\\nsaturation\" at 0.1,100000 center");
        gnuplot.plot("pc-Sw");

        gnuplot.resetAll();

        const auto [krw, krn] = sampleRelPerms(makeFluidMatrixInteraction(pcKrSwCurve_), sw); // test wrapped law
        Gnuplot::addRelPerms(gnuplot, sw, krw, krn, "kr-Sw", "w lp");
        gnuplot.plot("kr");
    }

    /*!
     * \brief Applies the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param globalPos The global position
     */
    Scalar permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Defines the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return finePorosity_;
        else
            return coarsePorosity_;
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
    { return FluidSystem::H2OIdx; }

private:
    bool isFineMaterial_(const GlobalPosition &globalPos) const
    { return globalPos[dimWorld-1] > layerBottom_; }

    Scalar fineK_;
    Scalar coarseK_;
    Scalar layerBottom_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    const PcKrSwCurve pcKrSwCurve_;

    bool plotFluidMatrixInteractions_;
};

} // end namespace Dumux

#endif
