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
 * \brief Test problem for the sequential 2p models
 */

#ifndef TEST_2P_SPATIALPARAMETERS_HH
#define TEST_2P_SPATIALPARAMETERS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>


namespace Dumux
{

template<class TypeTag>
class Test2PSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(Test2PSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(Test2PSpatialParams, SpatialParams, Test2PSpatialParams<TypeTag>);

// Set the material law
SET_PROP(Test2PSpatialParams, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedBrooksCorey<Scalar> RawMaterialLaw;
public:
    typedef EffToAbsLaw<RawMaterialLaw> type;
};
}

/** \todo Please doc me! */

template<class TypeTag>
class Test2PSpatialParams: public FVSpatialParams<TypeTag>
{
    typedef FVSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;

    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };
    typedef typename Grid::Traits::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    const FieldMatrix& intrinsicPermeabilityAtPos(const GlobalPosition& globalPos) const
    {
        if (isLensOne(globalPos))
            return permLenses_;
        else if (isLensTwo(globalPos))
            return permLenses_;
        else if (isLensThree(globalPos))
            return permLenses_;
        else
            return permBackground_;
    }

    Scalar porosity(const Element& element) const
    {
#if PROBLEM == 0
        return 0.2;
#elif PROBLEM == 1
        return 0.3;
#else
        return 0.4;
#endif
    }

    // return the brooks-corey context depending on the position
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        if (isLensOne(globalPos))
            return materialLawParamsLenses_;
        else if (isLensTwo(globalPos))
            return materialLawParamsLenses_;
        else if (isLensThree(globalPos))
            return materialLawParamsLenses_;
        else
            return materialLawParamsBackground_;
    }

    Test2PSpatialParams(const GridView& gridView) :
            ParentType(gridView), permBackground_(0), permLenses_(0),
            lensOneLowerLeft_(0), lensOneUpperRight_(0), lensTwoLowerLeft_(0), lensTwoUpperRight_(0), lensThreeLowerLeft_(0), lensThreeUpperRight_(0)
    {
#if PROBLEM == 0
        // residual saturations
        materialLawParamsBackground_.setSwr(0.2);
        materialLawParamsBackground_.setSnr(0.2);

        materialLawParamsLenses_.setSwr(0.2);
        materialLawParamsLenses_.setSnr(0.2);

        //parameters for Brooks-Corey law

        // entry pressures function
        materialLawParamsBackground_.setPe(0.);
        materialLawParamsLenses_.setPe(0.);
#elif PROBLEM == 1
        // residual saturations
        materialLawParamsBackground_.setSwr(0.);
        materialLawParamsBackground_.setSnr(0.);

        materialLawParamsLenses_.setSwr(0.);
        materialLawParamsLenses_.setSnr(0.);

        //parameters for Brooks-Corey law

        // entry pressures function
        materialLawParamsBackground_.setPe(5000.);
        materialLawParamsLenses_.setPe(5000.);
#else
        // residual saturations
        materialLawParamsBackground_.setSwr(0.);
        materialLawParamsBackground_.setSnr(0.);

        materialLawParamsLenses_.setSwr(0.);
        materialLawParamsLenses_.setSnr(0.);

        //parameters for Brooks-Corey law

        // entry pressures function
        materialLawParamsBackground_.setPe(0.);
        if (ParameterTree::tree().hasKey("SpatialParams.BackgroundEntryPressure"))
        {
            materialLawParamsBackground_.setPe(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, BackgroundEntryPressure));
        }

        materialLawParamsLenses_.setPe(0.);
        if (ParameterTree::tree().hasKey("SpatialParams.LenseEntryPressure"))
        {
            materialLawParamsLenses_.setPe(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LenseEntryPressure));
        }
#endif

        // Brooks-Corey shape parameters

#if PROBLEM == 2
        materialLawParamsBackground_.setLambda(3);
        if (ParameterTree::tree().hasKey("SpatialParams.BackgroundLambda"))
        {
            materialLawParamsBackground_.setLambda(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, BackgroundLambda));
        }
        materialLawParamsLenses_.setLambda(2);
        if (ParameterTree::tree().hasKey("SpatialParams.LenseLambda"))
        {
            materialLawParamsLenses_.setLambda(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LenseLambda));
        }
#else
        materialLawParamsBackground_.setLambda(2);
        materialLawParamsLenses_.setLambda(2);
#endif

#if PROBLEM == 0
        permBackground_[0][0] = 1e-7;
        permBackground_[1][1] = 1e-7;
        permLenses_[0][0] = 1e-7;
        permLenses_[1][1] = 1e-7;
#else
        permBackground_[0][0] = 1e-10;
        permBackground_[1][1] = 1e-10;
        permLenses_[0][0] = 1e-10;
        permLenses_[1][1] = 1e-10;
#endif


#if PROBLEM == 2
        if (ParameterTree::tree().hasKey("SpatialParams.BackgroundPermeabilityXX"))
        {
            permBackground_[0][0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, BackgroundPermeabilityXX);
        }
        if (ParameterTree::tree().hasKey("SpatialParams.BackgroundPermeabilityXY"))
        {
            permBackground_[0][1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, BackgroundPermeabilityXY);
        }
        if (ParameterTree::tree().hasKey("SpatialParams.BackgroundPermeabilityYX"))
        {
            permBackground_[1][0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, BackgroundPermeabilityYX);
        }
        if (ParameterTree::tree().hasKey("SpatialParams.BackgroundPermeabilityYY"))
        {
            permBackground_[1][1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, BackgroundPermeabilityYY);
        }

        if (ParameterTree::tree().hasKey("SpatialParams.LensPermeabilityXX"))
        {
            permLenses_[0][0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LensPermeabilityXX);
        }
        if (ParameterTree::tree().hasKey("SpatialParams.LensPermeabilityXY"))
        {
            permLenses_[0][1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LensPermeabilityXY);
        }
        if (ParameterTree::tree().hasKey("SpatialParams.LensPermeabilityYX"))
        {
            permLenses_[1][0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LensPermeabilityYX);
        }
        if (ParameterTree::tree().hasKey("SpatialParams.LensPermeabilityYY"))
        {
            permLenses_[1][1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LensPermeabilityYY);
        }

        if (ParameterTree::tree().hasKey("SpatialParams.LensOneLowerLeftX"))
        {
            lensOneLowerLeft_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LensOneLowerLeftX);
        }

        if (ParameterTree::tree().hasKey("SpatialParams.LensOneUpperRightX"))
        {
            lensOneUpperRight_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LensOneUpperRightX);
        }

        if (ParameterTree::tree().hasKey("SpatialParams.LensTwoLowerLeftX"))
        {
            lensTwoLowerLeft_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LensTwoLowerLeftX);
        }
        if (ParameterTree::tree().hasKey("SpatialParams.LensTwoUpperRightX"))
        {
            lensTwoUpperRight_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LensTwoUpperRightX);
        }

        if (ParameterTree::tree().hasKey("SpatialParams.LensThreeLowerLeftX"))
        {
            lensThreeLowerLeft_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LensThreeLowerLeftX);
        }
        if (ParameterTree::tree().hasKey("SpatialParams.LensThreeUpperRightX"))
        {
            lensThreeUpperRight_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LensThreeUpperRightX);
        }

        if (ParameterTree::tree().hasKey("SpatialParams.LensOneLowerLeftY"))
        {
            lensOneLowerLeft_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LensOneLowerLeftY);
        }
        if (ParameterTree::tree().hasKey("SpatialParams.LensOneUpperRightY"))
        {
            lensOneUpperRight_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LensOneUpperRightY);
        }

        if (ParameterTree::tree().hasKey("SpatialParams.LensTwoLowerLeftY"))
        {
            lensTwoLowerLeft_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LensTwoLowerLeftY);
        }
        if (ParameterTree::tree().hasKey("SpatialParams.LensTwoUpperRightY"))
        {
            lensTwoUpperRight_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LensTwoUpperRightY);
        }

        if (ParameterTree::tree().hasKey("SpatialParams.LensThreeLowerLeftY"))
        {
            lensThreeLowerLeft_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LensThreeLowerLeftY);
        }
        if (ParameterTree::tree().hasKey("SpatialParams.LensThreeUpperRightY"))
        {
            lensThreeUpperRight_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LensThreeUpperRightY);
        }
#endif
    }

private:

    bool isLensOne(const GlobalPosition& globalPos) const
    {
        for (int i = 0; i < dim; i++)
        {
            if (globalPos[i] < lensOneLowerLeft_[i] || globalPos[i] > lensOneUpperRight_[i])
            {
                return false;
            }
        }
        return true;
    }
    bool isLensTwo(const GlobalPosition& globalPos) const
    {
        for (int i = 0; i < dim; i++)
        {
            if (globalPos[i] < lensTwoLowerLeft_[i] || globalPos[i] > lensTwoUpperRight_[i])
            {
                return false;
            }
        }
        return true;
    }
    bool isLensThree(const GlobalPosition& globalPos) const
    {
        for (int i = 0; i < dim; i++)
        {
            if (globalPos[i] < lensThreeLowerLeft_[i] || globalPos[i] > lensThreeUpperRight_[i])
            {
                return false;
            }
        }
        return true;
    }

    MaterialLawParams materialLawParamsBackground_;
    MaterialLawParams materialLawParamsLenses_;
    FieldMatrix permBackground_;
    FieldMatrix permLenses_;
    GlobalPosition lensOneLowerLeft_;
    GlobalPosition lensOneUpperRight_;
    GlobalPosition lensTwoLowerLeft_;
    GlobalPosition lensTwoUpperRight_;
    GlobalPosition lensThreeLowerLeft_;
    GlobalPosition lensThreeUpperRight_;
};

} // end namespace
#endif
