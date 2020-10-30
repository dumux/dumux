/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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

#include <dumux/common/properties.hh>
#include <dumux/material/spatialparams/sequentialfv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

namespace Dumux
{

template<class TypeTag>
class Test2PSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
namespace TTag {
struct Test2PSpatialParams {};
}

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Test2PSpatialParams> { using type = Test2PSpatialParams<TypeTag>; };

}

// forward declaration
template<class Scalar>
class LinearMaterialDefault;
class LinearMaterial;

/*!
 * \ingroup SequentialTwoPTests
 * \brief Test problem for the sequential 2p models
 */
template<class TypeTag>
class Test2PSpatialParams: public SequentialFVSpatialParams<TypeTag>
{
    using ParentType = SequentialFVSpatialParams<TypeTag>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CoordScalar = typename Grid::ctype;

    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };
    using Element = typename Grid::Traits::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FieldMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

public:

    static constexpr bool pcSwCurveIsLinear()
    {
        return std::is_same_v<PcKrSwCurve, LinearMaterial> || std::is_same_v<PcKrSwCurve, LinearMaterialDefault>;
    }

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

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     *
     * \param globalPos The global coordinates for the given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        if (isLensOne(globalPos) || isLensTwo(globalPos) || isLensThree(globalPos))
            return makeFluidMatrixInteraction(*pcKrSwCurveLenses_);
        else
            return makeFluidMatrixInteraction(*pcKrSwCurveBackground_);
    }

    Test2PSpatialParams(const Problem& problem) :
            ParentType(problem), permBackground_(0), permLenses_(0),
            lensOneLowerLeft_(0), lensOneUpperRight_(0), lensTwoLowerLeft_(0), lensTwoUpperRight_(0), lensThreeLowerLeft_(0), lensThreeUpperRight_(0)
    {
#if PROBLEM == 0
        typename PcKrSwCurve::BasicParams params(0/*pe*/, 2/*lambda*/);
        typename PcKrSwCurve::EffToAbsParams effToAbsParams(0.2/*swr*/, 0.2/*snr*/);
        pcKrSwCurveBackground_ = std::make_unique<PcKrSwCurve>(params, effToAbsParams);
        pcKrSwCurveLenses_ = std::make_unique<PcKrSwCurve>(params, effToAbsParams);

#elif PROBLEM == 1
        typename PcKrSwCurve::BasicParams params(5000/*pe*/, 2/*lambda*/);
        typename PcKrSwCurve::EffToAbsParams effToAbsParams(0/*swr*/, 0/*snr*/);
        pcKrSwCurveBackground_ = std::make_unique<PcKrSwCurve>(params, effToAbsParams);
        pcKrSwCurveLenses_ = std::make_unique<PcKrSwCurve>(params, effToAbsParams);

#else
        typename PcKrSwCurve::BasicParams paramsBackground(getParam<Scalar>("SpatialParams.BackgroundEntryPressure", 0.0),
                                                           getParam<Scalar>("SpatialParams.BackgroundLambda", 3.0));
        typename PcKrSwCurve::BasicParams paramsLenses(getParam<Scalar>("SpatialParams.LenseEntryPressure", 0.0),
                                                           getParam<Scalar>("SpatialParams.LenseLambda", 2.0));
        typename PcKrSwCurve::EffToAbsParams effToAbsParams(0/*swr*/, 0/*snr*/);
        pcKrSwCurveBackground_ = std::make_unique<PcKrSwCurve>(paramsBackground, effToAbsParams);
        pcKrSwCurveLenses_ = std::make_unique<PcKrSwCurve>(paramsLenses, effToAbsParams);
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

    std::unique_ptr<PcKrSwCurve> pcKrSwCurveBackground_;
    std::unique_ptr<PcKrSwCurve> pcKrSwCurveLenses_;

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
