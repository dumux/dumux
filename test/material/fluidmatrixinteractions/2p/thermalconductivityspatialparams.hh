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
 * \brief Definition of the spatial parameters for the thermal conductivity tests.
 */
#ifndef DUMUX_THERMAL_CONDUCTIVITY_SPATIAL_PARAMS_HH
#define DUMUX_THERMAL_CONDUCTIVITY_SPATIAL_PARAMS_HH

#include <dumux/io/plotthermalconductivitymodel.hh>
#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/2p2c/implicit/model.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class ThermalConductivitySpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(ThermalConductivitySpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(ThermalConductivitySpatialParams, SpatialParams, ThermalConductivitySpatialParams<TypeTag>);

// Set the material law parameterized by absolute saturations
SET_TYPE_PROP(ThermalConductivitySpatialParams,
              MaterialLaw,
              EffToAbsLaw<RegularizedBrooksCorey<typename GET_PROP_TYPE(TypeTag, Scalar)> >);

// Define whether to open a gnuplot window
NEW_PROP_TAG(OutputOpenPlotWindow);
SET_BOOL_PROP(ThermalConductivitySpatialParams, OutputOpenPlotWindow, false);
}

/*!
 * \ingroup MaterialTestProblems
 * \brief Definition of the spatial parameters for the thermal conductivity tests.
 */
template<class TypeTag>
class ThermalConductivitySpatialParams : public ImplicitSpatialParams<TypeTag>
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
    ThermalConductivitySpatialParams(const GridView &gridView)
        : ParentType(gridView)
    {
        porosity_ = 0.3;
        rhoSolid_ = 2700.0;
        lambdaSolid_ = 2.8;

        // residual saturations
        materialParams_.setSwr(0.2);
        materialParams_.setSnr(0.0);

        // parameters for the Brooks-Corey law
        materialParams_.setPe(1e4);
        materialParams_.setLambda(2.0);
    }

    /*!
     * \brief This is called from the problem and creates a gnuplot output
     *        of e.g the pc-Sw curve
     */
    void plotMaterialLaw()
    {
        GnuplotInterface<Scalar> gnuplot;
        gnuplot.setOpenPlotWindow(GET_PARAM_FROM_GROUP(TypeTag, bool, Output, OpenPlotWindow));
        PlotThermalConductivityModel<TypeTag> plotThermalConductivityModel_(293.15, 1e5);
        std::string fileName = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Conductivity, File);
        plotThermalConductivityModel_.addlambdaeffcurve(gnuplot, porosity_, rhoSolid_, lambdaSolid_,
                                                        0.0, 1.0, fileName);
        gnuplot.plot("lambda_eff");
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
    { return 1e-10; }

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
    { return porosity_; }


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
    { return materialParams_; }

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
    { return 790; /* specific heat capacity of granite [J / (kg K)] */ }

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
    { return rhoSolid_; /* density of granite [kg/m^3] */ }

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
    { return lambdaSolid_; /* [W/(m K) */ }

private:
    Scalar porosity_;
    Scalar lambdaSolid_;
    Scalar rhoSolid_;

    MaterialLawParams materialParams_;
};

}

#endif
