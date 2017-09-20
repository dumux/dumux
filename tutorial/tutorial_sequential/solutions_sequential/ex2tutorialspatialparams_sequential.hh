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
 * \brief spatial parameters for the sequential tutorial
 */
#ifndef DUMUX_EX2TUTORIAL_SPATIAL_PARAMS_SEQUENTIAL_HH
#define DUMUX_EX2TUTORIAL_SPATIAL_PARAMS_SEQUENTIAL_HH


#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class Ex2TutorialSpatialParamsSequential;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(Ex2TutorialSpatialParamsSequential);

// Set the spatial parameters
SET_TYPE_PROP(Ex2TutorialSpatialParamsSequential, SpatialParams,
        Ex2TutorialSpatialParamsSequential<TypeTag>); /*@\label{tutorial-sequential:set-spatialparameters}@*/

// Set the material law
SET_PROP(Ex2TutorialSpatialParamsSequential, MaterialLaw)
{
private:
    // material law typedefs
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedBrooksCorey<Scalar> RawMaterialLaw;
public:
    typedef EffToAbsLaw<RawMaterialLaw> type;
};
}

//! Definition of the spatial parameters for the sequential tutorial

template<class TypeTag>
class Ex2TutorialSpatialParamsSequential: public FVSpatialParams<TypeTag>
{
    typedef FVSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;

    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld};

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    /**
     * \brief Intrinsic permeability tensor K \f$[m^2]\f$ depending
     *        on the position in the domain
     *
     *  \param globalPos The global position in the domain.
     *
     *  Alternatively, the function intrinsicPermeabilityAtPos(const GlobalPosition& globalPos) could be
     *  defined, where globalPos is the vector including the global coordinates of the finite volume.
     */
    const FieldMatrix& intrinsicPermeabilityAtPos(const GlobalPosition& globalPos) const
    {
        if(((globalPos[0]>25)&&(globalPos[0]<75))&&((globalPos[1]>15)&&(globalPos[1]<35)))
            return KLense_;
        else
            return K_;
    }

    /**
     * \brief Define the porosity \f$[-]\f$ of the porous medium depending
     *        on the position in the domain
     *
     *  \param globalPos The global position in the domain.
     *
     *  Alternatively, the function porosityAtPos(const GlobalPosition& globalPos) could be
     *  defined, where globalPos is the vector including the global coordinates of the finite volume.
     */
    double porosityAtPos(const GlobalPosition& globalPos) const
    {
        if(((globalPos[0]>25)&&(globalPos[0]<75))&&((globalPos[1]>15)&&(globalPos[1]<35)))
        return 0.15;
        else
        return 0.2;
    }

    /*! \brief Return the parameter object for the material law (i.e. Brooks-Corey)
     *         depending on the position in the domain
     *
     *  \param globalPos The global position in the domain.
     *
     *  Alternatively, the function materialLawParamsAtPos(const GlobalPosition& globalPos)
     *  could be defined, where globalPos is the vector including the global coordinates of
     *  the finite volume.
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        if(((globalPos[0]>25)&&(globalPos[0]<75))&&((globalPos[1]>15)&&(globalPos[1]<35)))
        return materialParamsLense_;
        else
        return materialParams_;
    }

    //! Constructor
    Ex2TutorialSpatialParamsSequential(const GridView& gridView)
    : ParentType(gridView), K_(0), KLense_(0)
    {
        for (int i = 0; i < dim; i++)
                K_[i][i] = 1e-7;
        for (int i = 0; i < dim; i++)
                KLense_[i][i] = 1e-9;


        //set residual saturations
            materialParams_.setSwr(0.0);                /*@\label{tutorial-coupled:setLawParams}@*/
            materialParams_.setSnr(0.0);
            materialParamsLense_.setSwr(0.0);                /*@\label{tutorial-coupled:setLawParams}@*/
            materialParamsLense_.setSnr(0.0);

            //parameters of Brooks & Corey Law
            materialParams_.setPe(1000.0);
            materialParams_.setLambda(1.8);
            materialParamsLense_.setPe(1500.0);
            materialParamsLense_.setLambda(2.0);
    }

private:
    FieldMatrix K_;
    FieldMatrix KLense_;

    MaterialLawParams materialParams_;                 /*@\label{tutorial-coupled:matParamsObject}@*/
    MaterialLawParams materialParamsLense_;
   };

} // end namespace
#endif
