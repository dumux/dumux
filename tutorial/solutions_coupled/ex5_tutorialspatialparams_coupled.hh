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
#ifndef DUMUX_TUTORIAL_EX5_SPATIAL_PARAMS_COUPLED_HH
#define DUMUX_TUTORIAL_EX5_SPATIAL_PARAMS_COUPLED_HH

// include parent spatialparameters
#include <dumux/material/spatialparams/implicit.hh>

// include material laws
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh> /*@\label{tutorial-coupled:rawLawInclude}@*/
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>

namespace Dumux {
//forward declaration
template<class TypeTag>
class Ex5TutorialSpatialParamsCoupled;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(Ex5TutorialSpatialParamsCoupled);/*@\label{tutorial-coupled:define-spatialparameters-typetag}@*/

// Set the spatial parameters
SET_TYPE_PROP(Ex5TutorialSpatialParamsCoupled, SpatialParams,
        Dumux::Ex5TutorialSpatialParamsCoupled<TypeTag>); /*@\label{tutorial-coupled:set-spatialparameters}@*/

// Set the material law
SET_PROP(Ex5TutorialSpatialParamsCoupled, MaterialLaw)
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
class Ex5TutorialSpatialParamsCoupled: public ImplicitSpatialParams<TypeTag> /*@\label{tutorial-coupled:tutorialSpatialParameters}@*/
{
    // Get informations for current implementation via property system
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum
    {
        dim = Grid::dimension
    };

    // Get object types for function arguments
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;

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
    { return K_; }

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
    { return 0.2; }

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
        return materialParams_;
    }

    // constructor
    Ex5TutorialSpatialParamsCoupled(const GridView& gridView) :
        ImplicitSpatialParams<TypeTag>(gridView),
        K_(0)
    {
        //set main diagonal entries of the permeability tensor to a value
        //setting to one value means: isotropic, homogeneous
        for (int i = 0; i < dim; i++)
            K_[i][i] = 1e-7;

        //set residual saturations
        materialParams_.setSwr(0.0);                /*@\label{tutorial-coupled:setLawParams}@*/
        materialParams_.setSnr(0.0);

        //parameters of Brooks & Corey Law
        materialParams_.setPe(500.0);
        materialParams_.setLambda(2);
    }

private:
    Dune::FieldMatrix<Scalar, dim, dim> K_;
    // Object that holds the values/parameters of the selected material law.
    MaterialLawParams materialParams_;                 /*@\label{tutorial-coupled:matParamsObject}@*/
};
} // end namespace
#endif
