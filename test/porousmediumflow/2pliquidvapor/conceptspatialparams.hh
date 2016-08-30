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
 * \brief Definition of the spatial parameters for various steam injection problems
 */



#ifndef DUMUX_INJECTION_SPATIAL_PARAMS_HH
#define DUMUX_INJECTION_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class InjectionProblemSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(InjectionProblemSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(InjectionProblemSpatialParams, SpatialParams, Dumux::InjectionProblemSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(InjectionProblemSpatialParams, MaterialLaw)
{
 private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedVanGenuchten<Scalar> EffMaterialLaw;
 public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffMaterialLaw> type;
};

// This file is used for multiple problem setups which might require different input parameters.
// Therefore, default values are set here to assure functionality:

NEW_PROP_TAG(SpatialParamsWell);
NEW_PROP_TAG(SpatialParamsPermeability);

SET_STRING_PROP(InjectionProblemSpatialParams, SpatialParamsWell, "default");
SET_SCALAR_PROP(InjectionProblemSpatialParams, SpatialParamsPermeability, 1e-11);

}

/*!
 * \ingroup TwoPOneCNIModel
 *
 * \brief Definition of the spatial parameters for various steam injection problems
 */
template<class TypeTag>
class InjectionProblemSpatialParams : public ImplicitSpatialParams<TypeTag>
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

//     typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
//     enum {
//         wPhaseIdx = Indices::wPhaseIdx,
//         gPhaseIdx = Indices::gPhaseIdx
//     };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar,dimWorld> DimVector;


    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;


public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    InjectionProblemSpatialParams(const GridView &gridView)
        : ParentType(gridView)
    {

        eps_ = 1e-6;

        // intrinsic permeability
        permeability_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Permeability); //1e-11;//1e-16;

        // porosity
        porosity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Porosity);

        // heat conductivity of granite
        lambdaSolid_ = 2.8;

        // specific heat capacities
        heatCap_ = 850.;

        // set Van Genuchten Parameters
        materialParams_.setSwr(0.1);
        materialParams_.setSnr(0.0);
        materialParams_.setVgAlpha(0.0028);
        materialParams_.setVgn(2.0);

        // Anisotropy factor for K
        anisotropyFactor_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.KHorToKVert);

        // the vertical permeability distribution is specified by the name of a certain well  (Durlach - Ochs, 2006)
        well_ = GET_PARAM_FROM_GROUP(TypeTag,
                                  std::string,
                                  SpatialParams,
                                  Well);
    }

    ~InjectionProblemSpatialParams()
    {}



    /*!
     * \brief Update the spatial parameters with the flow solution
     *        after a timestep.
     *
     * \param globalSolution The global solution vector
     */
    void update(const SolutionVector &globalSolution)
    {
    };

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const FieldMatrix intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvElemGeom,
                                       int scvIdx) const
    {
        const GlobalPosition &pos = fvElemGeom.subContVol[scvIdx].global;


        FieldMatrix permMatrix;
        permMatrix = 0.0;
        Scalar ratio = anisotropyFactor_ ;

        permMatrix[0][0] =1;
        permMatrix[1][1] =1;
        permMatrix[2][2] = 1/ratio ;

   //distribution of K for well Br38
   if(well_ == "Br38"){
        if(pos[2] < 0.5 + eps_){
            permMatrix *= 1e-9;
            return permMatrix;
        }

        else if(pos[2] > 0.5 && pos[2] < 6 +eps_ ){
            permMatrix *= 5.5e-11;
            return permMatrix;
        }


        else {
            permMatrix *= 6.5e-12;
            return permMatrix;
        }
        }

      //distribution of K for well I6
    else if(well_ == "I6"){
         if(pos[2] <= 1.5 /*+ eps_*/) {
             permMatrix *= 8e-10;
             return permMatrix;
         }

        else if(pos[2] > 1.5 && pos[2] <=2.5 /*+eps_*/ ) {
            permMatrix *= 9.3e-11;
             return permMatrix;
        }

       else if(pos[2] > 2.5 && pos[2] <= 3 && sqrt(pos[0]*pos[0] + pos[1]*pos[1]) > 1.5/*+eps_*/ ) {
           permMatrix *= 9.3e-11;
           return permMatrix;
       }

       else if(pos[2] > 2.5 && pos[2] <= 3 && sqrt(pos[0]*pos[0] + pos[1]*pos[1]) <= 1.5/*+eps_*/ ) {
           permMatrix *= 8e-11;
           return permMatrix;
       }

       else if(pos[2] > 3 && pos[2] < 5 /*+eps_ */) {
           permMatrix *= 2.7e-11;
            return permMatrix;
       }

       else if(pos[2] >= 5 && pos[2] <= 6 /*+eps_*/ ) {
           permMatrix *= 8e-12;
            return permMatrix;
       }

        else {
            permMatrix *= 5.3e-12;
            return  permMatrix;
        }
        }


   //distribution of K for well EK2
   else if(well_ == "EK2"){
        if(pos[2] < 2 + eps_) {
            permMatrix *= 3e-10;
            return permMatrix;
        }

        else if(pos[2] > 2 && pos[2] < 3 +eps_ ) {
            permMatrix *= 4e-11;
            return permMatrix;
        }

       else if(pos[2] > 3 && pos[2] < 4 +eps_ ) {
           permMatrix *= 1e-11;
            return permMatrix;
       }

       else if(pos[2] > 4 && pos[2] < 6 +eps_ ) {
           permMatrix *= 8.5e-12;
            return permMatrix;
       }

        else {
            permMatrix *= 5e-12;
            return permMatrix;
        }
        }
    else {
        permMatrix *= permeability_;
        return permMatrix; //default value
    }

    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    double porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    const int scvIdx) const
    {
        return porosity_;
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
        return materialParams_;
    }

    /*!
     * \brief Returns the heat capacity \f$[J/m^3 K]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the heat capacity needs to be defined
     */
    double heatCapacity(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const int scvIdx) const
    {
            return heatCap_ * 2650 // density of sand [kg/m^3]
                * (1 - porosity(element, fvGeometry, scvIdx));
    }



     Scalar solidHeatCapacity(const Element &element,
                             const FVElementGeometry &fvGeometry,
                             const int scvIdx) const
    {
        return heatCap_; // specific heat capacity of granite [J / (kg K)]
    }


    Scalar solidDensity(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const int scvIdx) const
    {
        return 2650; // density of granite [kg/m^3]
    }


     Scalar solidThermalConductivity(const Element &element,
                                    const FVElementGeometry &fvGeometry,
                                    const int scvIdx) const
        {
          return lambdaSolid_;
        }

private:

    Scalar lambdaSolid_;
    Scalar eps_;
    std::string well_;
    Scalar permeability_;
    Scalar porosity_;
    Scalar heatCap_;
    Scalar anisotropyFactor_ ;
    MaterialLawParams materialParams_;
};

}

#endif
