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
class TwoPAntonioSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TwoPAntonioSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TwoPAntonioSpatialParams, SpatialParams, Dumux::TwoPAntonioSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(TwoPAntonioSpatialParams, MaterialLaw)
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
class TwoPAntonioSpatialParams : public ImplicitSpatialParams<TypeTag>
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


    TwoPAntonioSpatialParams(const GridView &gridView)
    : ParentType(gridView)
    {
        // episode index
        episode_ = 0;
        // intrinsic permeabilities [m^2]
        Kinit_ = Scalar(0.0); // init permeability
        KFault_ = Scalar(0.0); // Fault permeability
        KMatrix_ = Scalar(0.0); // Matrix permeability
        KCaprock_ = Scalar(0.0); // Matrix permeability
        for (int i = 0; i < dim; i++){
            Kinit_[i][i] = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.Kinit);
            KFault_[i][i] = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.KFault);
            KMatrix_[i][i] = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.KMatrix);
            KCaprock_[i][i] = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.KCaprock);
        }

        // Tolerance for fault zone geometry
        spatialTolarance_ = GET_RUNTIME_PARAM(TypeTag, Scalar,FailureParameters.SpatialTolerance);

        // porosities [-]
        phiMatrix_  = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.PhiMatrix);
        phiCaprock_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.PhiCaprock);

        // rock density [kg/m^3]
        rockDensity_ = 2260.0;

        // Young's modulus of the pure elastic model (= only a spring)
        Ematrix_ = GET_RUNTIME_PARAM(TypeTag, Scalar,ElasticParameters.Ematrix);
        Efault_ = GET_RUNTIME_PARAM(TypeTag, Scalar,ElasticParameters.Efault);

        // Poisson's ratio [-]
        nuMatrix_ = GET_RUNTIME_PARAM(TypeTag, Scalar,ElasticParameters.nuMatrix);
        nuFault_ = GET_RUNTIME_PARAM(TypeTag, Scalar,ElasticParameters.nuFault);

        // Lame parameters [Pa]
        lambdaMatrix_ = (Ematrix_ * nuMatrix_) / ((1 + nuMatrix_)*(1 - 2 * nuMatrix_));
        lambdaFault_  = (Efault_ * nuFault_)  / ((1 + nuFault_ )*(1 - 2 * nuFault_));

        muMatrix_ = Ematrix_ / (2 * (1 + nuMatrix_));
        muFault_  = Efault_  / (2 * (1 + nuFault_));

        poreCompFault_ = GET_RUNTIME_PARAM(TypeTag, Scalar,ElasticParameters.PoreCompFault);
        poreCompMatrix_ = GET_RUNTIME_PARAM(TypeTag, Scalar,ElasticParameters.PoreCompMatrix);
        poreCompCaprock_ = GET_RUNTIME_PARAM(TypeTag, Scalar,ElasticParameters.PoreCompCaprock);

        faultAngle_ = 90 + GET_RUNTIME_PARAM(TypeTag, Scalar,FailureParameters.FaultAngle);

        frictionAngleFault_ = GET_RUNTIME_PARAM(TypeTag, Scalar,FailureParameters.FrictionAngleFault) * M_PI / 180;
        cohesionFault_ = GET_RUNTIME_PARAM(TypeTag, Scalar,FailureParameters.CohesionFault);

        frictionAngleMatrix_ = GET_RUNTIME_PARAM(TypeTag, Scalar,FailureParameters.FrictionAngleMatrix) * M_PI / 180;
        cohesionMatrix_ = GET_RUNTIME_PARAM(TypeTag, Scalar,FailureParameters.CohesionMatrix);

        // Viscosity
        // viscosity_= 1.0285e15;//2.0187e11;
        deltaT_ = GET_RUNTIME_PARAM(TypeTag, Scalar,FailureParameters.DeltaT); // stress drop is defined with respect to a certain timelength deltaT [s], which is set here
        deltaSigma_=  GET_RUNTIME_PARAM(TypeTag, Scalar,FailureParameters.DeltaSigma);

        // given Van Genuchten m
        m_ = 0.457;

        n_ = 1.0 / (1.0 - m_);

        alpha_ = 1.0 / (19.9e3);

        // residual saturations
        MaterialParams_.setSwr(0.20);
        MaterialParams_.setSnr(0.05);

        MaterialParams_.setVgAlpha(alpha_);
        MaterialParams_.setVgn(n_);

        MaterialParams_.setKrwHighSw(0.99);

        //////////////////////////////////////

        MaterialParams2_.setSwr(0.30);
        MaterialParams2_.setSnr(0.05);

        MaterialParams2_.setVgAlpha(alpha_);
        MaterialParams2_.setVgn(n_);

        MaterialParams2_.setKrwHighSw(0.99);

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

    ~TwoPAntonioSpatialParams()
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

        if (isInFault_(globalPos, xFaultBottom_, yFaultBottom_, xFaultTop_, yFaultTop_))
        {
            return phiMatrix_;
        }

        else if (isCaprock_(globalPos,  xCaprockBottom_,  yCaprockBottom_, xCaprockTop_, yCaprockTop_))
        {
            return phiCaprock_;
        }
        else
            return phiMatrix_;

    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the soil
     *
     * \param globalPos The global position of the vertex
     */
    double porosityAtPos(const GlobalPosition& globalPos) const
    {
        if (isInFault_(globalPos, xFaultBottom_, yFaultBottom_, xFaultTop_, yFaultTop_))
        {
                return phiMatrix_;
        }

        else if (isCaprock_(globalPos,  xCaprockBottom_,  yCaprockBottom_, xCaprockTop_, yCaprockTop_))
        {
            return phiCaprock_;
        }
        else
            return phiMatrix_;

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
        if(episode_ < 1)
            return Kinit_; // intrinsic permeability applied during initialization
        else
        {
            if (isInFault_(globalPos, xFaultBottom_, yFaultBottom_, xFaultTop_, yFaultTop_))
            {
                return KFault_;
            }

            else if (isCaprock_(globalPos,  xCaprockBottom_,  yCaprockBottom_, xCaprockTop_, yCaprockTop_))
            {
                return KCaprock_;
            }
            else
                return KMatrix_;
        }
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

        if (isInFault_(globalPos, xFaultBottom_, yFaultBottom_, xFaultTop_, yFaultTop_))
        {
//             std::cout << "I am here!" << std::endl;
            param[0] = lambdaFault_;
            param[1] = muFault_;
            param[2] = poreCompFault_;
        }
//         else if (isInFault_(center, xFaultBottom2_, yFaultBottom2_, xFaultTop2_, yFaultTop2_))
//         {
// //             std::cout << "I am here!" << std::endl;
//             param[0] = lambdaFault_ * 0.5;
//             param[1] = muFault_ * 0.5;
//             param[2] = poreCompFault_ * 0.5;
//
//         }
//             else if (xCoord > 487.5 && xCoord < 512.5 &&
//                      yCoord > 987.5   && yCoord < 1012.5)
//             {
//                 param[0] = Efault_;
//                 param[1] = Bfault_;
//             }
        else if (isCaprock_(globalPos,  xCaprockBottom_,  yCaprockBottom_, xCaprockTop_, yCaprockTop_))
        {
            param[0] = lambdaMatrix_;
            param[1] = muMatrix_;
            param[2] = poreCompCaprock_;
        }

        else
        {
            param[0] = lambdaMatrix_;
            param[1] = muMatrix_;
            param[2] = poreCompMatrix_;
        }
//         }

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

        if (isInFault_(globalPos, xFaultBottom_, yFaultBottom_, xFaultTop_, yFaultTop_))
        {

            failureCurveParams[0] = frictionAngleFault_;
            failureCurveParams[1] = cohesionFault_;
        }
//             else if (xCoord > 487.5 && xCoord < 512.5 &&
// //                      yCoord < 50 &&
// //                      yCoord > 40 &&
//                      yCoord > 987.5   && yCoord < 1012.5)
//             {
//                 failureCurveParams[0] = frictionAngleFault_;
//                 failureCurveParams[1] = cohesionFault_;
//             }
        else
            failureCurveParams[0] = frictionAngleMatrix_;
            failureCurveParams[1] = cohesionMatrix_;
//         }

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

    const MaterialLawParams& materialLawParams2(const Element &element,
                                               const FVElementGeometry &fvGeometry,
                                               int scvIdx) const
    {
        return MaterialParams2_;
    }

private:
    bool isInFault_(const GlobalPosition &globalPos,
                    const Scalar xFaultBottom,
                    const Scalar yFaultBottom,
                    const Scalar xFaultTop,
                    const Scalar yFaultTop) const
    {
        // x-Values of cell centers at the bottom: 320.3, 322.8, 325.3, 327.8
        // x-Values of cell centers at the top:    672.2, 674.7, 677.2, 679.7

        Scalar xCoord = globalPos[0];
        Scalar yCoord = globalPos[1];

//         Scalar xFaultTop =  674.7;
//         Scalar yFaultTop =  1000.0;
//
//         Scalar xFaultBottom =  322.8;
//         Scalar yFaultBottom = -1000.0;

        Scalar tanAngle = (yFaultTop - yFaultBottom) / (xFaultTop - xFaultBottom);

        Scalar deltaXBottom = xCoord - xFaultBottom;
        Scalar deltaYBottom = yCoord - yFaultBottom;

        Scalar deltaXTop = xFaultTop - xCoord;
        Scalar deltaYTop = yFaultTop - yCoord;

        if (deltaYBottom < tanAngle * deltaXBottom + spatialTolarance_ &&
            deltaYBottom > tanAngle * deltaXBottom - spatialTolarance_ &&
            std::abs(yCoord) <  500 + eps_ &&
            deltaYBottom > yFaultBottom &&
            deltaYBottom < yFaultTop)
        {
            return true;
        }
        else if (deltaYTop < tanAngle * deltaXTop + spatialTolarance_ &&
            deltaYTop > tanAngle * deltaXTop - spatialTolarance_ &&
            std::abs(yCoord) <  500 + eps_ &&
            deltaYTop > yFaultBottom &&
            deltaYTop < yFaultTop)
        {
            return true;
        }
        else
            return false;
    }

    bool isCaprock_(const GlobalPosition &globalPos,
                    const Scalar xCaprockBottom,
                    const Scalar yCaprockBottom,
                    const Scalar xCaprockTop,
                    const Scalar yCaprockTop) const
    {
        // x-Values of cell centers at the bottom: 320.3, 322.8, 325.3, 327.8
        // x-Values of cell centers at the top:    672.2, 674.7, 677.2, 679.7

        Scalar xCoord = globalPos[0];
        Scalar yCoord = std::abs(globalPos[1]);

        if (xCoord < xCaprockTop + eps_ &&
            xCoord > xCaprockBottom - eps_  &&
            yCoord < yCaprockTop + eps_ &&
            yCoord > yCaprockBottom - eps_  )
        {
            return true;
        }
        else
            return false;
    }



    Scalar spatialTolarance_;
    Dune::FieldMatrix<Scalar,dim,dim> KFault_, KMatrix_, Kinit_, KCaprock_;
    Scalar layerBottom_;
    Scalar rockDensity_;
    Scalar phiMatrix_, phiCaprock_;
    Scalar lambda_;
    Scalar mu_;
    Scalar Ematrix_, Efault_;
    Scalar nuMatrix_, nuFault_;
    Scalar lambdaMatrix_ ,lambdaFault_;
    Scalar muMatrix_, muFault_;
    Scalar poreCompFault_, poreCompMatrix_, poreCompCaprock_;
    Scalar faultAngle_;
    Scalar frictionAngleMatrix_, frictionAngleFault_;
    Scalar cohesionMatrix_, cohesionFault_;
    //Scalar viscosity_;
    Scalar deltaT_;
    Scalar deltaSigma_;
    Scalar BrooksCoreyLambda_, m_, alpha_, n_;
    MaterialLawParams MaterialParams_;
    MaterialLawParams MaterialParams2_;
    static constexpr Scalar eps_ = 3e-6;
    int episode_;
    bool anyFailure_;

    bool plotFluidMatrixInteractions_;

    // x-Values of cell centers at the bottom: 320.3, 322.8, 325.3, 327.8
    // x-Values of cell centers at the top:    672.2, 674.7, 677.2, 679.7
    Scalar xFaultTop_ =  674.7;
    Scalar yFaultTop_ =  1000.0;

    Scalar xFaultTop2_ =  677.2;
    Scalar yFaultTop2_ =  1000.0;

    Scalar xFaultBottom_ =  322.8;
    Scalar yFaultBottom_ = -1000.0;

    Scalar xFaultBottom2_ =  325.3;
    Scalar yFaultBottom2_ = -1000.0;

    Scalar xCaprockTop_ =  2000;
    Scalar yCaprockTop_ =  200;

    Scalar xCaprockBottom_ =  0.0;
    Scalar yCaprockBottom_ =   50;

};
}
#endif
