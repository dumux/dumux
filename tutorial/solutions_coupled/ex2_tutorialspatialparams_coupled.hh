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
 * \brief The spatial parameters for the fully coupled tutorial problem
 *        which uses the twophase box model.
 */
#ifndef DUMUX_EX2_TUTORIAL_SPATIAL_PARAMS_COUPLED_HH
#define DUMUX_EX2_TUTORIAL_SPATIAL_PARAMS_COUPLED_HH

// include parent spatialparameters
#include <dumux/material/spatialparams/implicit.hh>

// include material laws
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh> /*@\label{tutorial-coupled:rawLawInclude}@*/
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>

namespace Dumux {
//forward declaration
template<class TypeTag>
class Ex2TutorialSpatialParamsCoupled;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(Ex2TutorialSpatialParamsCoupled);/*@\label{tutorial-coupled:define-spatialparameters-typetag}@*/

// Set the spatial parameters
SET_TYPE_PROP(Ex2TutorialSpatialParamsCoupled, SpatialParams,
        Dumux::Ex2TutorialSpatialParamsCoupled<TypeTag>); /*@\label{tutorial-coupled:set-spatialparameters}@*/

// Set the material law
SET_PROP(Ex2TutorialSpatialParamsCoupled, MaterialLaw)
{
private:
    // material law typedefs
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    // select material law to be used
    typedef RegularizedBrooksCorey<Scalar> RawMaterialLaw;     /*@\label{tutorial-coupled:rawlaw}@*/
public:
    // adapter for absolute law
    typedef EffToAbsLaw<RawMaterialLaw> type;   /*@\label{tutorial-coupled:eff2abs}@*/
};
}

/*!
 * \ingroup TwoPBoxModel
 *
 * \brief The spatial parameters for the fully coupled tutorial problem
 *        which uses the twophase box model.
 */
template<class TypeTag>
class Ex2TutorialSpatialParamsCoupled: public ImplicitSpatialParams<TypeTag> /*@\label{tutorial-coupled:tutorialSpatialParameters}@*/
{
    // Get informations for current implementation via property system
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum
    {
        dim = Grid::dimension,
    dimWorld = Grid::dimensionworld
    };

    // Get object types for function arguments
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    // get material law from property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    // determine appropriate parameters depending on selected materialLaw
    typedef typename MaterialLaw::Params MaterialLawParams;    /*@\label{tutorial-coupled:matLawObjectType}@*/

    /*! Intrinsic permeability tensor K \f$[m^2]\f$ depending
     *  on the position in the domain
     *
     *  \param element The finite volume element
     *  \param fvGeometry The finite-volume geometry in the box scheme
     *  \param scvIdx The local vertex index
     *
     *  Alternatively, the function intrinsicPermeabilityAtPos(const GlobalPosition& globalPos)
     *  could be defined, where globalPos is the vector including the global coordinates
     *  of the finite volume.
     */
    const Dune::FieldMatrix<Scalar, dim, dim> &intrinsicPermeability(const Element &element, /*@\label{tutorial-coupled:permeability}@*/
                                                    const FVElementGeometry &fvGeometry,
                                                    const int scvIdx) const
    {
        GlobalPosition globalPos = element.geometry().corner(scvIdx);

        if (globalPos[0] > 25 && globalPos[0] < 75 && globalPos[1] > 15 && globalPos[1] < 35)
            return KLense_;
        else
            return K_;
     }

    /*! Defines the porosity \f$[-]\f$ of the porous medium depending
     * on the position in the domain
     *
     *  \param element The finite volume element
     *  \param fvGeometry The finite-volume geometry in the box scheme
     *  \param scvIdx The local vertex index
     *
     *  Alternatively, the function porosityAtPos(const GlobalPosition& globalPos)
     *  could be defined, where globalPos is the vector including the global coordinates
     *  of the finite volume.
     */
    Scalar porosity(const Element &element,                    /*@\label{tutorial-coupled:porosity}@*/
            const FVElementGeometry &fvGeometry,
            const int scvIdx) const
    {
        GlobalPosition globalPos = element.geometry().corner(scvIdx);

        if (globalPos[0] > 25 && globalPos[0] < 75 && globalPos[1] > 15 && globalPos[1] < 35)
            return porosityLense_;
        else
            return porosity_;
    }

    /*! Returns the parameter object for the material law (i.e. Brooks-Corey)
     *  depending on the position in the domain
     *
     *  \param element The finite volume element
     *  \param fvGeometry The finite-volume geometry in the box scheme
     *  \param scvIdx The local vertex index
     *
     *  Alternatively, the function materialLawParamsAtPos(const GlobalPosition& globalPos)
     *  could be defined, where globalPos is the vector including the global coordinates
     *  of the finite volume.
     */
    const MaterialLawParams& materialLawParams(const Element &element,            /*@\label{tutorial-coupled:matLawParams}@*/
                                               const FVElementGeometry &fvGeometry,
                                               const int scvIdx) const
    {
        GlobalPosition globalPos = element.geometry().corner(scvIdx);

        if (globalPos[0] > 25 && globalPos[0] < 75 && globalPos[1] > 15 && globalPos[1] < 35)
            return materialParamsLense_;
        else
            return materialParams_;
    }

    // constructor
    Ex2TutorialSpatialParamsCoupled(const GridView& gridView) :
        ImplicitSpatialParams<TypeTag>(gridView),
        K_(0), KLense_(0)
    {
        //set main diagonal entries of the permeability tensor to a value
        //setting to one value means: isotropic, homogeneous
        for (int i = 0; i < dim; i++)
        {
            K_[i][i] = 1e-7;
            KLense_[i][i] = 1e-9;
        }

        porosity_ = 0.2;
        porosityLense_ = 0.15;

        //set residual saturations
        materialParams_.setSwr(0.0);                /*@\label{tutorial-coupled:setLawParams}@*/
        materialParams_.setSnr(0.0);

        //parameters of Brooks & Corey Law
        materialParams_.setPe(1000.0);
        materialParams_.setLambda(1.8);

        //set residual saturations for lense
        materialParamsLense_.setSwr(0.0);                /*@\label{tutorial-coupled:setLawParams}@*/
        materialParamsLense_.setSnr(0.0);

        //parameters of Brooks & Corey Law
        materialParamsLense_.setPe(1500.0);
        materialParamsLense_.setLambda(2.0);
    }

private:
    Dune::FieldMatrix<Scalar, dim, dim> K_;
    Dune::FieldMatrix<Scalar, dim, dim> KLense_;
    Scalar porosity_;
    Scalar porosityLense_;
    // Object that holds the values/parameters of the selected material law.
    MaterialLawParams materialParams_;                 /*@\label{tutorial-coupled:matParamsObject}@*/
    MaterialLawParams materialParamsLense_;
};
} // end namespace
#endif
