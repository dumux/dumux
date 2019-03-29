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
 * \brief The spatial parameters for the ElRichards_TestProblem which uses the
 *        linear elastic two-phase model
 */
#ifndef DUMUX_ELRICHARDSPARAMETERS_HH
#define DUMUX_ELRICHARDSPARAMETERS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>//!!!!!here there is no equivalent folder for richards, therefore kept them. it should be ok because richards is a special case for 2p
//#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
//#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include "regularizedvangenuchten.hh"
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/geomechanics/elrichards/model.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class ElRichardsSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(ElRichardsSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(ElRichardsSpatialParams, SpatialParams, ElRichardsSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(ElRichardsSpatialParams, MaterialLaw)
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
 * \ingroup ElRichardsBoxModel
 * \brief The spatial parameters for the ElRichards_TestProblem which uses the
 *        linear elastic Richards model
 */
template<class TypeTag>
class ElRichardsSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld,
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> DimMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;

public:
    //get the material law from the property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;


    ElRichardsSpatialParams(const GridView &gridView)
    : ParentType(gridView)
    {
        // episode index
        episode_ = 0;
        // intrinsic permeabilities [m^2]
        Kinit_ = Scalar(0.0); // init permeability
        K_ = Scalar(0.0); // permeability
        k_ = Scalar(0.0); // conductivity.khodam
        ks_ = GET_RUNTIME_PARAM(TypeTag, Scalar, HydraulicParameters.ks);//khodam [m/s] we assume this is for standard temp and presure


        for (int i = 0; i < dim; i++){ //we consider a value almost KInit=100*K to make it larger for faster initialization
        Kinit_ = 1.E-12; //[m²]//this larger value is just to speedup the initializaion simulation
        //K_[i][i] = 1.E-14; //[m²]//this is the permeability which depends on the soil properties. later on in localresidual.hh for flux calculation,  effective permeability will be calculated from this K*krw which krw depends on saturation condition and was defind in vangenuchten.hh
        K_ = (ks_*1.E-3)/(1000*9.81);//khodam [m/s] dynamic viscosity and density are defined from simpleh2o.hh in dumux/material/component , and didnot make them pressure and temp dependent, becasue Ks is constant therefore I consider the standard tem,press for it. k=ks*Krw
//         krw_ =  MaterialLaw::krw(matParams,);
//         k_[i][i] = ks_*krw_;

        }

        // porosities [-]
        phi_ = 0.46;
//         Scalar phi_ = GET_RUNTIME_PARAM(TypeTag, Scalar, HydraulicParameters.phi);

        // rock density [kg/m^3]
//      rockDensity_ = 2000;

        rockDensity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, FailureParameters.rockDensity);
        // Young's modulus [Pa]
      //E_ = 1e4;
        E_ = GET_RUNTIME_PARAM(TypeTag, Scalar, FailureParameters.E);




        // Poisson's ratio [-]
        //nu_ = 0.33;
        nu_ = GET_RUNTIME_PARAM(TypeTag, Scalar, FailureParameters.nu);
        // Lame parameters [Pa]
        lambda_ = (E_ * nu_) / ((1 + nu_)*(1 - 2 * nu_)); //elasticity
        mu_ = E_ / (2 * (1 + nu_));


        // given Van Genuchten m
 //   m_ = 1.37;

      Scalar n_ = GET_RUNTIME_PARAM(TypeTag, Scalar, HydraulicParameters.n);//has been calculated in vangenuchtenparams.hh
      Scalar l_ = GET_RUNTIME_PARAM(TypeTag, Scalar, HydraulicParameters.l);
      Scalar alpha_ = GET_RUNTIME_PARAM(TypeTag, Scalar, HydraulicParameters.alpha);
      Scalar Swr = GET_RUNTIME_PARAM(TypeTag, Scalar, HydraulicParameters.residualWaterContent);

     //Friction Angle & Cohesion has been used in model.hh
          theta_ = GET_RUNTIME_PARAM(TypeTag, Scalar, FailureParameters.FrictionAngle)*M_PI / 6;
          s0_ = GET_RUNTIME_PARAM(TypeTag, Scalar, FailureParameters.Cohesion);


        // Brooks Corey lambda
        //BrooksCoreyLambda_ = m_ / (1 - m_) * (1 - pow(0.5,1/m_));

        // residual saturations
        MaterialParams_.setSwr(Swr);
        //MaterialParams_.setSnr(0.034);

        // parameters for the van Genuchten law
        MaterialParams_.setVgn(n_);
        MaterialParams_.setVgl(l_);
        MaterialParams_.setVgAlpha(alpha_);
//         MaterialParams_.setKrw(krw_); //khodam


     }

    ~ElRichardsSpatialParams()
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
    const Scalar intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvGeometry,
                                       int scvIdx) const
    {
        if(episode_ <= 1)
           // return Kinit_; // intrinsic permeability applied during initialization
            return Kinit_;
        else
             return K_; // intrinsic permeability
    }


    /*!
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
     * \brief Define the friction angle \f$[-]\f$ of the soil
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    Scalar theta(const Element &element,
                    const FVElementGeometry &fvGeometry) const
    {
        const GlobalPosition& globalPos = element.geometry().center();

//         if (isTopLayer_(globalPos))
//             return thetaOver_;
//         else
//              return thetaUnder_;
        return theta_;
    }

    /*!
     * \brief Define the friction angle \f$[-]\f$ of the soil
     *
     * \param globalPos The global position of the vertex
     */
    Scalar theta(const GlobalPosition& globalPos) const
    {
//         if (isTopLayer_(globalPos))
//               return thetaOver_;
//         return thetaUnder_;
        return theta_;
    }

    /*!
     * \brief Define the cohesion \f$[-]\f$ of the soil
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    Scalar S0(const Element &element,
                    const FVElementGeometry &fvGeometry) const
    {
        const GlobalPosition& globalPos = element.geometry().center();

//         if (isTopLayer_(globalPos))
//             return s0Over_;
//         return s0Under_;
        return s0_;
    }

    /*!
     * \brief Define the cohesion \f$[-]\f$ of the soil
     *
     * \param globalPos The global position of the vertex
     */
    Scalar S0(const GlobalPosition& globalPos) const
    {
//         if (isTopLayer_(globalPos))
//             return s0Over_;
//         return s0Under_;
        return s0_;
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

        param[0] = lambda_; //elasticity
        param[1] = mu_;

        return param;
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
    Scalar K_, Kinit_;
    Scalar k_, ks_;//khodam
    Scalar rockDensity_;
    Scalar phi_, phiInit_;
    Scalar lambda_;
    Scalar mu_;
    Scalar E_;
    Scalar nu_;
//     Scalar l_;
//     Scalar n_;
//     Scalar alpha_;
    Scalar krw_;
    //Scalar BrooksCoreyLambda_;
    Scalar m_;
    Scalar theta_, s0_;
    MaterialLawParams MaterialParams_;
    static constexpr Scalar eps_ = 3e-6;
    int episode_;


};
}
#endif
