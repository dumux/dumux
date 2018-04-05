// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPTwoCTests
 * \brief Definition of the spatial parameters for the water-air problem.
 */
//forward declaration
template<class TypeTag>
class WaterAirSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(WaterAirSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(WaterAirSpatialParams, SpatialParams, WaterAirSpatialParams<TypeTag>);


// Set the material law parameterized by absolute saturations
SET_TYPE_PROP(WaterAirSpatialParams,
              MaterialLaw,
              EffToAbsLaw<RegularizedBrooksCorey<typename GET_PROP_TYPE(TypeTag, Scalar)> >);
}

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief Definition of the spatial parameters for the water-air problem
 */
template<class TypeTag>
class WaterAirSpatialParams : public FVSpatialParams<TypeTag>
{
    using ParentType = FVSpatialParams<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using CoordScalar = typename GridView::ctype;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

public:
    //! export the type used for the permeability
    using PermeabilityType = Scalar;
    //! export the type used for the material law
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    WaterAirSpatialParams(const Problem& problem) : ParentType(problem)
    {
        layerBottom_ = 22.0;

        // intrinsic permeabilities
        fineK_ = 1e-13;
        coarseK_ = 1e-12;

        // porosities
        finePorosity_ = 0.3;
        coarsePorosity_ = 0.3;

        // heat conductivity of granite
        lambdaSolid_ = 2.8;

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

        gnuplot.resetAll();
        using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);
        PlotEffectiveDiffusivityModel<Scalar, EffDiffModel> plotEffectiveDiffusivityModel;
        plotEffectiveDiffusivityModel.adddeffcurve(gnuplot, finePorosity_, 0.0, 1.0, "fine");
        plotEffectiveDiffusivityModel.adddeffcurve(gnuplot, coarsePorosity_, 0.0, 1.0, "coarse");
        gnuplot.plot("deff");

        gnuplot.resetAll();
        using ThermCondModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);
        using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
        PlotThermalConductivityModel<Scalar, ThermCondModel, FluidSystem> plotThermalConductivityModel;
        plotThermalConductivityModel.addlambdaeffcurve(gnuplot, finePorosity_, 2700.0, lambdaSolid_, 0.0, 1.0, "fine");
        plotThermalConductivityModel.addlambdaeffcurve(gnuplot, coarsePorosity_, 2700.0, lambdaSolid_, 0.0, 1.0, "coarse");
        gnuplot.plot("lambdaeff");
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
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global positio
     */
    Scalar solidHeatCapacityAtPos(const GlobalPosition& globalPos) const
    { return 790; /*specific heat capacity of granite [J / (kg K)]*/ }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidDensityAtPos(const GlobalPosition& globalPos) const
    { return 2700; /*density of granite [kg/m^3]*/ }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * \param globalPos The global position
     */
    Scalar solidThermalConductivityAtPos(const GlobalPosition& globalPos) const
    { return lambdaSolid_; }

private:
    bool isFineMaterial_(const GlobalPosition &globalPos) const
    { return globalPos[dimWorld-1] > layerBottom_; }

    Scalar fineK_;
    Scalar coarseK_;
    Scalar layerBottom_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    // heat conductivity of the solid material only
    Scalar lambdaSolid_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    bool plotFluidMatrixInteractions_;
};

} // end namespace Dumux

#endif
