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

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/ploteffectivediffusivitymodel.hh>
#include <dumux/io/plotmateriallaw.hh>
#include <dumux/io/plotthermalconductivitymodel.hh>
#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTwoCModel
 * \brief Definition of the spatial parameters for the water-air problem
 */
template<class FVGridGeometry, class Scalar>
class WaterAirSpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar,
                         WaterAirSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar,
                                       WaterAirSpatialParams<FVGridGeometry, Scalar>>;

    static constexpr int dimWorld = GridView::dimensionworld;

    using EffectiveLaw = RegularizedBrooksCorey<Scalar>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! export the type used for the permeability
    using PermeabilityType = Scalar;
    //! export the type used for the material law
    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;

    /*!
     * \brief The constructor
     *
     * \param fvGridGeometry The finite volume grid geometry
     */
    WaterAirSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry) : ParentType(fvGridGeometry)
    {
        layerBottom_ = 22.0;

        // intrinsic permeabilities
        fineK_ = 1e-13;
        coarseK_ = 1e-12;

        // porosities
        finePorosity_ = 0.3;
        coarsePorosity_ = 0.3;

        // residual saturations
        fineMaterialParams_.setSwr(0.2);
        fineMaterialParams_.setSnr(0.0);
        coarseMaterialParams_.setSwr(0.2);
        coarseMaterialParams_.setSnr(0.0);

        // parameters for the Brooks-Corey law
        fineMaterialParams_.setPe(1e4);
        coarseMaterialParams_.setPe(1e4);
        fineMaterialParams_.setLambda(2.0);
        coarseMaterialParams_.setLambda(2.0);

        plotFluidMatrixInteractions_ = getParam<bool>("Output.PlotFluidMatrixInteractions");
    }

    /*!
     * \brief This is called from the problem and creates a gnuplot output
     *        of e.g the pc-Sw curve
     */
    void plotMaterialLaw()
    {
        PlotMaterialLaw<Scalar, MaterialLaw> plotMaterialLaw;
        GnuplotInterface<Scalar> gnuplot(plotFluidMatrixInteractions_);
        gnuplot.setOpenPlotWindow(plotFluidMatrixInteractions_);
        plotMaterialLaw.addpcswcurve(gnuplot, fineMaterialParams_, 0.2, 1.0, "fine", "w lp");
        plotMaterialLaw.addpcswcurve(gnuplot, coarseMaterialParams_, 0.2, 1.0, "coarse", "w l");
        gnuplot.setOption("set xrange [0:1]");
        gnuplot.setOption("set label \"residual\\nsaturation\" at 0.1,100000 center");
        gnuplot.plot("pc-Sw");

        gnuplot.resetAll();
        plotMaterialLaw.addkrcurves(gnuplot, fineMaterialParams_, 0.2, 1.0, "fine");
        plotMaterialLaw.addkrcurves(gnuplot, coarseMaterialParams_, 0.2, 1.0, "coarse");
        gnuplot.plot("kr");
    }

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
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
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
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
     * \brief return the parameter object for the Brooks-Corey material law which depends on the position
     *
     * \param globalPos The global position
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return fineMaterialParams_;
        else
            return coarseMaterialParams_;
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The position of the center of the element
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

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    bool plotFluidMatrixInteractions_;
};

} // end namespace Dumux

#endif
