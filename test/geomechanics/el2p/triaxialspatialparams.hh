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
#ifndef DUMUX_TRIAXIALELTWOPSPARAMETERS_HH
#define DUMUX_TRIAXIALTWOPSPARAMETERS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/geomechanics/el2p/model.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TriaxialEl2PSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TriaxialEl2PSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TriaxialEl2PSpatialParams, SpatialParams, Dumux::TriaxialEl2PSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(TriaxialEl2PSpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedBrooksCorey<Scalar> EffectiveLaw;
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
class TriaxialEl2PSpatialParams : public ImplicitSpatialParams<TypeTag>
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

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> DimMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;

public:
    //get the material law from the property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;


    TriaxialEl2PSpatialParams(const GridView &gridView)
    : ParentType(gridView)
    {
        // episode index
        episode_ = 0;
        // intrinsic permeabilities [m^2]
        Kinit_ = Scalar(0.0); // init permeability
        KFault_ = Scalar(0.0); // Fault permeability
        KMatrix_ = Scalar(0.0); // Matrix permeability
        for (int i = 0; i < dim; i++){
            Kinit_[i][i] = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.Kinit);
            KFault_[i][i] = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.KFault);
            KMatrix_[i][i] = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.KMatrix);
        }

        // Tolerance for fault zone geometry
        spatialTolarance_ = GET_RUNTIME_PARAM(TypeTag, Scalar,FailureParameters.SpatialTolerance);

        // porosities [-]
        phi_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.Phi);

        // rock density [kg/m^3]
        rockDensity_ = 2700.0;

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
        // Brooks Corey lambda
        BrooksCoreyLambda_ = m_ / (1 - m_) * (1 - pow(0.5,1/m_));

        // residual saturations
        MaterialParams_.setSwr(0.3);
        MaterialParams_.setSnr(0.05);

        // parameters for the Brooks Corey law
        MaterialParams_.setPe(1.99e4);
        MaterialParams_.setLambda(BrooksCoreyLambda_);


     }

    ~TriaxialEl2PSpatialParams()
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

    void setAnyFailure(const bool& anyFailure)
    {
        anyFailure_ = anyFailure;
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
        if(episode_ <= 1)
            return Kinit_; // intrinsic permeability applied during initialization
        else
        {
            GlobalPosition center = element.geometry().center();

            if (isInFault_(center))
            {
                return KFault_;
            }
//             else if (xCoord > 487.5 && xCoord < 512.5 &&
//                      yCoord > 987.5   && yCoord < 1012.5)
//             {
//                 return KFault_;
//             }
            else
                return KMatrix_;
        }
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
            return phi_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the soil
     *
     * \param globalPos The global position of the vertex
     */
    double porosity(const GlobalPosition& globalPos) const
    {
            return phi_;
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
    const Dune::FieldVector<Scalar,2> lameParams(const Element &element,
                                           const FVElementGeometry &fvGeometry,
                                           int scvIdx) const
    {
        // Lame parameters
        Dune::FieldVector<Scalar, 2> param;

        GlobalPosition center = element.geometry().center();

        if (isInFault_(center))
        {
//             std::cout << "I am here!" << std::endl;
            param[0] = lambdaFault_;
            param[1] = muFault_;
        }
//             else if (xCoord > 487.5 && xCoord < 512.5 &&
//                      yCoord > 987.5   && yCoord < 1012.5)
//             {
//                 param[0] = Efault_;
//                 param[1] = Bfault_;
//             }
        else
        {
            param[0] = lambdaMatrix_;
            param[1] = muMatrix_;
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

        GlobalPosition center = element.geometry().center();

        if (isInFault_(center))
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

private:
    bool isInFault_(const GlobalPosition &globalPos) const
    {
        Scalar xCoord = globalPos[0];
        Scalar yCoord = globalPos[1];

        Scalar xCenterFault = 0.5;
        Scalar yCenterFault = 0.5;

//         Scalar xFaultDistanceCenterToBothTips = 87.5; // x-distance in m to both tips of the fault
        Scalar xFaultDistanceCenterToBothTips = 100; // x-distance in m to both tips of the fault

        Scalar yFaultBottom = yCenterFault - tan(55 * M_PI / 180) * xFaultDistanceCenterToBothTips;
        Scalar yFaultTop = yCenterFault + tan(55 * M_PI / 180) * xFaultDistanceCenterToBothTips;

        Scalar deltaX = xCenterFault - xCoord;

        if(
            (std::abs(deltaX) < xFaultDistanceCenterToBothTips + eps_ &&
            yCoord < yCenterFault - tan(55 * M_PI / 180) * deltaX + spatialTolarance_ &&
            yCoord > yCenterFault - tan(55 * M_PI / 180) * deltaX - spatialTolarance_ &&
//                 yCoord < 50 &&
//                 yCoord > 40 &&
            yCoord > yFaultBottom &&
            yCoord < yFaultTop)
            // for Heterogeneous
            &&
            (!(xCoord < 0.21 && xCoord > 0.09 && yCoord > 0.00 && yCoord < 0.11))
            &&
            (!(xCoord < 0.91 && xCoord > 0.79 && yCoord > 0.89 && yCoord < 1.01)) )
            return true;
//         else if (xCoord < 50 &&
//             xCoord > 45 &&
//             yCoord > 45 &&
//             yCoord < 50)
//             return true;
        else
            return false;
    }

    Scalar spatialTolarance_;
    Dune::FieldMatrix<Scalar,dim,dim> KFault_, KMatrix_, Kinit_;
    Scalar layerBottom_;
    Scalar rockDensity_;
    Scalar phi_, phiInit_;
    Scalar lambda_;
    Scalar mu_;
    Scalar Ematrix_, Efault_;
    Scalar nuMatrix_, nuFault_;
    Scalar lambdaMatrix_ ,lambdaFault_;
    Scalar muMatrix_, muFault_;
    Scalar faultAngle_;
    Scalar frictionAngleMatrix_, frictionAngleFault_;
    Scalar cohesionMatrix_, cohesionFault_;
    //Scalar viscosity_;
    Scalar deltaT_;
    Scalar deltaSigma_;
    Scalar BrooksCoreyLambda_, m_;
    MaterialLawParams MaterialParams_;
    static constexpr Scalar eps_ = 3e-6;
    int episode_;
    bool anyFailure_;

};
}
#endif
