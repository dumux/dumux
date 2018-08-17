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
 * \brief The spatial parameters for the El2P_TestProblem which uses the
 *        linear elastic two-phase model
 */
#ifndef DUMUX_TWOPANTONIOSPARAMETERS_HH
#define DUMUX_TWOPANTONIOSPARAMETERS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/io/plotmateriallaw.hh>

// #include <dumux/geomechanics/viscoel2p/model.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TwoPSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TwoPSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TwoPSpatialParams, SpatialParams, Dumux::TwoPSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(TwoPSpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedVanGenuchten<Scalar> EffectiveLaw;
public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffectiveLaw> type;
};
}
/*!
 * \ingroup ElTwoPBoxModel
 * \brief The spatial parameters for the ViscoEl2P_TestProblem which uses the
 *        linear elastic two-phase model
 */
template<class TypeTag>
class TwoPSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld,
    };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };


    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> DimMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;

public:
    //get the material law from the property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;


    TwoPSpatialParams(const GridView &gridView)
    : ParentType(gridView)
    {
        // episode index
        episode_ = 0;
        // intrinsic permeabilities [m^2]
        Kinit_ = Scalar(0.0); // init permeability
//         K_ = Scalar(0,0); //Permeability
        k_ = Scalar(0.0); // conductivity.khodam
        ks_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.ks);//khodam [m/s] we assume this is for standard temp and presure

        for (int i = 0; i < dim; i++){
            Kinit_[i][i] = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.Kinit);
            K_ = (ks_*1.E-3)/(1000*9.81);//khodam [m/s] dynamic viscosity and density are defined from simpleh2o.hh in dumux/material/component , and didnot make them pressure and temp dependent, becasue Ks is constant therefore I consider the standard tem,press for it. k=ks*Krw

        }

        // Tolerance for fault zone geometry
        spatialTolarance_ = GET_RUNTIME_PARAM(TypeTag, Scalar,FailureParameters.SpatialTolerance);

        // porosities [-]
        phi_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.Phi);

        // rock density [kg/m^3]
        rockDensity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, FailureParameters.rockDensity);

        // Young's modulus of the pure elastic model (= only a spring)
        E_ = GET_RUNTIME_PARAM(TypeTag, Scalar,ElasticParameters.E);

        // Poisson's ratio [-]
        nu_ = GET_RUNTIME_PARAM(TypeTag, Scalar,ElasticParameters.nu);

        // Lame parameters [Pa]
        lambda_ = (E_ * nu_) / ((1 + nu_)*(1 - 2 * nu_));

        mu_ = E_ / (2 * (1 + nu_));

        poreComp_ = GET_RUNTIME_PARAM(TypeTag, Scalar,ElasticParameters.PoreComp);


//         faultAngle_ = 90 + GET_RUNTIME_PARAM(TypeTag, Scalar,FailureParameters.FaultAngle);

        frictionAngle_ = GET_RUNTIME_PARAM(TypeTag, Scalar,FailureParameters.FrictionAngle) * M_PI / 180;
        cohesion_ = GET_RUNTIME_PARAM(TypeTag, Scalar,FailureParameters.Cohesion);



        // Viscosity
        // viscosity_= 1.0285e15;//2.0187e11;
        deltaT_ = GET_RUNTIME_PARAM(TypeTag, Scalar,FailureParameters.DeltaT); // stress drop is defined with respect to a certain timelength deltaT [s], which is set here
        deltaSigma_=  GET_RUNTIME_PARAM(TypeTag, Scalar,FailureParameters.DeltaSigma);

        // given Van Genuchten m
        n_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.n);

        m_ = 1.0 - (1.0 / n_);

        alpha_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.alpha);
        Scalar swr_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.swr); //residual water content
        Scalar snr_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.snr);

        // residual saturations
        MaterialParams_.setSwr(swr_);
        MaterialParams_.setSnr(snr_);

        MaterialParams_.setVgAlpha(alpha_);
        MaterialParams_.setVgn(n_);

//         MaterialParams_.setKrwHighSw(0.99);


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
        plotMaterialLaw.addpcswcurve(gnuplot, MaterialParams_, 0.2, 1.0, "fine", "w lp");
        gnuplot.setOption("set xrange [0:1]");
        gnuplot.setOption("set label \"residual\\nsaturation\" at 0.1,100000 center");
        gnuplot.plot("pc-Sw");

//         gnuplot.resetAll();
//         plotMaterialLaw.addkrcurves(gnuplot, fineMaterialParams_, 0.2, 1.0, "fine");
//         plotMaterialLaw.addkrcurves(gnuplot, coarseMaterialParams_, 0.2, 1.0, "coarse");
//         gnuplot.plot("kr");
//
//         gnuplot.resetAll();
//         PlotEffectiveDiffusivityModel<TypeTag> plotEffectiveDiffusivityModel;
//         plotEffectiveDiffusivityModel.adddeffcurve(gnuplot, finePorosity_, 0.0, 1.0, "fine");
//         plotEffectiveDiffusivityModel.adddeffcurve(gnuplot, coarsePorosity_, 0.0, 1.0, "coarse");
//         gnuplot.plot("deff");
//
//         gnuplot.resetAll();
//         PlotThermalConductivityModel<TypeTag> plotThermalConductivityModel;
//         plotThermalConductivityModel.addlambdaeffcurve(gnuplot, finePorosity_, 2700.0, lambdaSolid_, 0.0, 1.0, "fine");
//         plotThermalConductivityModel.addlambdaeffcurve(gnuplot, coarsePorosity_, 2700.0, lambdaSolid_, 0.0, 1.0, "coarse");
//         gnuplot.plot("lambdaeff");
    }

    ~TwoPSpatialParams()
    {}

    /*!
     * \brief This function sets the private variable episode_ to the current episode index
     * which is checked in the hydraulic parameter functions to identify if we are still in the
     * initialization run (episode_ == 1)
     *
     * \param episode The episode index
     */
    void setEpisode(const int& episode)
    {
        episode_ = episode;
        std::cout<< "episode set to: "<< episode_<<std::endl;
    }

    /*!dumux efunktion
     * \brief Define the porosity \f$[-]\f$ of the soil
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    double porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    {
        GlobalPosition globalPos = element.geometry().center();

            return phi_;

    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the soil
     *
     * \param globalPos The global position of the vertex
     */
    double porosityAtPos(const GlobalPosition& globalPos) const
    {
            return phi_;
    }
    /*!
     * \brief Apply the intrinsic permeability tensor \f$[m^2]\f$ to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     *
     * During the initialization period the intrinsic permeability can be set to a larger
     * value in order to accelerate the calculation of the hydrostatic pressure distribution.
     */
    const DimMatrix intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvGeometry,
                                       int scvIdx) const
    {
        GlobalPosition globalPos = element.geometry().center();

        return intrinsicPermeabilityAtPos(globalPos);
    }

    const DimMatrix intrinsicPermeabilityAtPos(const GlobalPosition& globalPos) const
    {
        if(episode_ <= 1)
           // return Kinit_; // intrinsic permeability applied during initialization
            return Kinit_;
        else
             return K_; // intrinsic permeability
    }

    void setAnyFailure(const bool& anyFailure)
    {
        anyFailure_ = anyFailure;
    }

    /*!
     * \brief Define the density \f$[kg/m^3]\f$ of the rock
     *
     * \param element The finite element
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    const Scalar rockDensity(const Element &element,
                                        int scvIdx) const
    {
        return rockDensity_;
    }

    /*!
     * \brief Define the density \f$[kg/m^3]\f$ of the rock
     *
     * \param globalPos The global position of the vertex
     */
    const Scalar rockDensity(const GlobalPosition &globalPos) const
    {
        return rockDensity_;
    }

    /*!
     * \brief Define the Lame parameters \f$[Pa]\f$ linear elastic rock
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    const Dune::FieldVector<Scalar,3> lameParams(const Element &element,
                                           const FVElementGeometry &fvGeometry,
                                           int scvIdx) const
    {
        // Lame parameters
        Dune::FieldVector<Scalar, 3> param;

        GlobalPosition globalPos = element.geometry().center();


        param[0] = lambda_; //elasticity
        param[1] = mu_;

        return param;
    }

    /*!
     * \brief Define the parameters for failureCurve
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    const Dune::FieldVector<Scalar,2> failureCurveParams(const Element &element,
                                           const FVElementGeometry &fvGeometry,
                                           int scvIdx) const
    {
        // Lame parameters
        Dune::FieldVector<Scalar, 2> failureCurveParams;

        GlobalPosition globalPos = element.geometry().center();



            failureCurveParams[0] = frictionAngle_;
            failureCurveParams[1] = cohesion_;


        return failureCurveParams;
    }

    const Dune::FieldVector<Scalar,2> viscoParams(const Element &element,
                            const FVElementGeometry &fvGeometry,
                            int scvIdx) const
        {
        Dune::FieldVector<Scalar, 2> visco;
        //visco[0]=viscosity_;
        visco[0]= deltaT_;
        visco[1]= deltaSigma_;
        return visco;
        }


    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-Sw, pc-Sw, etc.).
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return the material parameters object
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const FVElementGeometry &fvGeometry,
                                               int scvIdx) const
    {
        return MaterialParams_;
    }

    const MaterialLawParams& materialLawParams(const GlobalPosition& globalPos) const//I added this then I can read the materialLawParams via spatialparams in problem to calculate sw just based on globalPos and no need for fvGeometry, element and scvIdx. then I dont need to read via model, as in localoperatopr.hh, and similar to porosirty I just read it via spatialparams.
    {
        return MaterialParams_;
    }

private:

    Scalar spatialTolarance_;
    Dune::FieldMatrix<Scalar,dim,dim>  Kinit_;
    Scalar K_;
    Scalar layerBottom_;
    Scalar rockDensity_;
    Scalar lambda_;
    Scalar mu_;
    Scalar poreComp_;
    Scalar frictionAngle_;
    Scalar cohesion_;
    //Scalar viscosity_;
    Scalar deltaT_;
    Scalar deltaSigma_;
    Scalar BrooksCoreyLambda_, m_, alpha_, n_;
    Scalar k_, ks_;//khodam
    Scalar phi_, phiInit_;
    Scalar E_;
    Scalar nu_;
    Scalar l_;


    MaterialLawParams MaterialParams_;
    MaterialLawParams MaterialParams2_;
    static constexpr Scalar eps_ = 3e-6;
    int episode_;
    bool anyFailure_;

    bool plotFluidMatrixInteractions_;



};
}
#endif
