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
 * \brief spatial parameters for the 2p 3d test
 */
#ifndef TEST_3D2P_SPATIALPARAMETERS_HH
#define TEST_3D2P_SPATIALPARAMETERS_HH

#include <dumux/material/spatialparams/sequentialfv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

namespace Dumux
{
/*!
 * \ingroup SequentialTwoPTests
 * \brief spatial parameters for the 2p 3d test
 */
//forward declaration
template<class TypeTag>
class Test3d2pSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
namespace TTag {
struct Test3d2pSpatialParams {};
}

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Test3d2pSpatialParams> { using type = Test3d2pSpatialParams<TypeTag>; };
}

/*!
 *
 * \ingroup IMPESTests
 * \brief spatial parameters for the 2p test using MPFAL 3D method
 */
template<class TypeTag>
class Test3d2pSpatialParams: public SequentialFVSpatialParams<TypeTag>
{
    using ParentType = SequentialFVSpatialParams<TypeTag>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CoordScalar = typename Grid::ctype;

    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld, numEq=1};
    using Element = typename Grid::Traits::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using LocalPosition = Dune::FieldVector<CoordScalar, dim>;
    using FieldMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

public:

    void update (Scalar saturationW, const Element& element)
    {
    }

    const FieldMatrix& intrinsicPermeabilityAtPos (const GlobalPosition& globalPos) const
    {
        return constPermeability_;
    }

    double porosity(const Element& element) const
    {
#if PROBLEM == 1
        return 0.3;
#else
        return 0.2;
#endif
    }


    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     *
     * \param globalPos The global coordinates for the given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        return makeFluidMatrixInteraction(*pcKrSwCurve_);
    }


    Test3d2pSpatialParams(const Problem& problem)
    : ParentType(problem), constPermeability_(0)
    {
        // parameters for the Brooks-Corey Law
#if PROBLEM == 1
        typename PcKrSwCurve::BasicParams params(5000/*pe*/, 2/*lambda*/);
        typename PcKrSwCurve::EffToAbsParams effToAbsParams(0.0/*swr*/, 0.0/*snr*/);
#else
        typename PcKrSwCurve::BasicParams params(0/*pe*/, 2/*lambda*/);
        typename PcKrSwCurve::EffToAbsParams effToAbsParams(0.2/*swr*/, 0.2/*snr*/);
#endif
        pcKrSwCurve_ = std::make_unique<PcKrSwCurve>(params, effToAbsParams);

#if PROBLEM == 1
        for(int i = 0; i < dim; i++)
        {
            constPermeability_[i][i] = 1e-10;
        }
#else
        for(int i = 0; i < dim; i++)
        {
            constPermeability_[i][i] = 1e-7;
        }
#endif
    }

private:
    FieldMatrix constPermeability_;
    std::unique_ptr<const PcKrSwCurve> pcKrSwCurve_;

};

} // end namespace
#endif
