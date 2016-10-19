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
 *
 * \brief Definition of the spatial parameters for the water-air problem.
 */
#ifndef DUMUX_WATER_AIR_SPATIAL_PARAMS_HH
#define DUMUX_WATER_AIR_SPATIAL_PARAMS_HH

#include <dumux/io/ploteffectivediffusivitymodel.hh>
#include <dumux/io/plotmateriallaw.hh>
#include <dumux/io/plotthermalconductivitymodel.hh>
#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/2p2c/implicit/model.hh>

namespace Dumux
{

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
class WaterAirSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    WaterAirSpatialParams(const GridView &gridView)
        : ParentType(gridView)
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

        plotFluidMatrixInteractions_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Output,
                                                                    PlotFluidMatrixInteractions);
    }

    /*!
     * \brief This is called from the problem and creates a gnuplot output
     *        of e.g the pc-Sw curve
     */
    void plotMaterialLaw()
    {
        PlotMaterialLaw<TypeTag> plotMaterialLaw;
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
        PlotEffectiveDiffusivityModel<TypeTag> plotEffectiveDiffusivityModel;
        plotEffectiveDiffusivityModel.adddeffcurve(gnuplot, finePorosity_, 0.0, 1.0, "fine");
        plotEffectiveDiffusivityModel.adddeffcurve(gnuplot, coarsePorosity_, 0.0, 1.0, "coarse");
        gnuplot.plot("deff");

        gnuplot.resetAll();
        PlotThermalConductivityModel<TypeTag> plotThermalConductivityModel;
        plotThermalConductivityModel.addlambdaeffcurve(gnuplot, finePorosity_, 2700.0, lambdaSolid_, 0.0, 1.0, "fine");
        plotThermalConductivityModel.addlambdaeffcurve(gnuplot, coarsePorosity_, 2700.0, lambdaSolid_, 0.0, 1.0, "coarse");
        gnuplot.plot("lambdaeff");
    }

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const Scalar intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvGeometry,
                                       const int scvIdx) const
    {
        const GlobalPosition &globalPos = fvGeometry.subContVol[scvIdx].global;
        if (isFineMaterial_(globalPos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    const int scvIdx) const
    {
        const GlobalPosition &globalPos = fvGeometry.subContVol[scvIdx].global;
        if (isFineMaterial_(globalPos))
            return finePorosity_;
        else
            return coarsePorosity_;
    }


    /*!
     * \brief return the parameter object for the Brooks-Corey material law which depends on the position
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const FVElementGeometry &fvGeometry,
                                               const int scvIdx) const
    {
        const GlobalPosition &globalPos = fvGeometry.subContVol[scvIdx].global;
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
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar solidHeatCapacity(const Element &element,
                             const FVElementGeometry &fvGeometry,
                             const int scvIdx) const
    {
        return 790; // specific heat capacity of granite [J / (kg K)]
    }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar solidDensity(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const int scvIdx) const
    {
        return 2700; // density of granite [kg/m^3]
    }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the heat capacity needs to be defined
     */
    Scalar solidThermalConductivity(const Element &element,
                                    const FVElementGeometry &fvGeometry,
                                    const int scvIdx) const
    {
        return lambdaSolid_;
    }

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

}

#endif
