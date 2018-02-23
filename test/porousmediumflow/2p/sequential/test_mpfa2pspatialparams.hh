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
 * \ingroup SequentialTwoPTests
 * \brief Test problem for the sequential 2p models
 */

#ifndef TEST_2P_SPATIALPARAMETERS_HH
#define TEST_2P_SPATIALPARAMETERS_HH

#include <dumux/material/spatialparams/sequentialfv.hh>
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
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using RawMaterialLaw = RegularizedBrooksCorey<Scalar>;
public:
    using type = EffToAbsLaw<RawMaterialLaw>;
};
}

/*!
 * \ingroup SequentialTwoPTests
 * \brief Test problem for the sequential 2p models \todo Please doc me!
 */
template<class TypeTag>
class Test2PSpatialParams: public SequentialFVSpatialParams<TypeTag>
{
    using ParentType = SequentialFVSpatialParams<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename Grid::ctype;

    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };
    using Element = typename Grid::Traits::template Codim<0>::Entity;

    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using FieldMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

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

    Test2PSpatialParams(const Problem& problem) :
            ParentType(problem), permBackground_(0), permLenses_(0),
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

        // entry pressures
        materialLawParamsBackground_.setPe(getParam<Scalar>("SpatialParams.BackgroundEntryPressure", 0.0));
        materialLawParamsLenses_.setPe(getParam<Scalar>("SpatialParams.LenseEntryPressure", 0.0));
#endif

        // Brooks-Corey shape parameters
#if PROBLEM == 2
        materialLawParamsBackground_.setLambda(getParam<Scalar>("SpatialParams.BackgroundLambda", 3.0));
        materialLawParamsLenses_.setLambda(getParam<Scalar>("SpatialParams.LenseLambda", 2.0));
#else
        materialLawParamsBackground_.setLambda(2.0);
        materialLawParamsLenses_.setLambda(2.0);
#endif

#if PROBLEM == 0
        permBackground_[0][0] = 1e-7;
        permBackground_[1][1] = 1e-7;
        permLenses_[0][0] = 1e-7;
        permLenses_[1][1] = 1e-7;
#elif PROBLEM == 1
        permBackground_[0][0] = 1e-10;
        permBackground_[1][1] = 1e-10;
        permLenses_[0][0] = 1e-10;
        permLenses_[1][1] = 1e-10;
#else
        permBackground_[0][0] = getParam<Scalar>("SpatialParams.BackgroundPermeabilityXX", 1e-10);
        permBackground_[0][1] = getParam<Scalar>("SpatialParams.BackgroundPermeabilityXY", 0.0);
        permBackground_[1][0] = getParam<Scalar>("SpatialParams.BackgroundPermeabilityYX", 0.0);
        permBackground_[1][1] = getParam<Scalar>("SpatialParams.BackgroundPermeabilityYY", 1e-10);
        permLenses_[0][0] = getParam<Scalar>("SpatialParams.LensPermeabilityXX", 1e-10);
        permLenses_[0][1] = getParam<Scalar>("SpatialParams.LensPermeabilityXY", 0.0);
        permLenses_[1][0] = getParam<Scalar>("SpatialParams.LensPermeabilityYX", 0.0);
        permLenses_[1][1] = getParam<Scalar>("SpatialParams.LensPermeabilityYY", 1e-10);

        lensOneLowerLeft_ = getParam<GlobalPosition>("SpatialParams.LensOneLowerLeft", GlobalPosition(0.0));
        lensOneUpperRight_ = getParam<GlobalPosition>("SpatialParams.LensOneUpperRight", GlobalPosition(0.0));
        lensTwoLowerLeft_ = getParam<GlobalPosition>("SpatialParams.LensTwoLowerLeft", GlobalPosition(0.0));
        lensTwoUpperRight_ = getParam<GlobalPosition>("SpatialParams.LensTwoUpperRight", GlobalPosition(0.0));
        lensThreeLowerLeft_ = getParam<GlobalPosition>("SpatialParams.LensThreeLowerLeft", GlobalPosition(0.0));
        lensThreeUpperRight_ = getParam<GlobalPosition>("SpatialParams.LensThreeUpperRight", GlobalPosition(0.0));
#endif
    }

private:

    bool isLensOne(const GlobalPosition& globalPos) const
    {
        for (int i = 0; i < dim; i++)
        {
            if (globalPos[i] < lensOneLowerLeft_[i] - eps_ || globalPos[i] > lensOneUpperRight_[i] + eps_)
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
            if (globalPos[i] < lensTwoLowerLeft_[i] - eps_ || globalPos[i] > lensTwoUpperRight_[i] + eps_)
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
            if (globalPos[i] < lensThreeLowerLeft_[i] - eps_ || globalPos[i] > lensThreeUpperRight_[i] + eps_)
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

    static constexpr Scalar eps_ = 1e-6;
};

} // end namespace
#endif
