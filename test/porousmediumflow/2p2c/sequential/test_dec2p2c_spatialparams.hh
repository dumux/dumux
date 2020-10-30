// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * \ingroup SequentialTwoPTwoCTests
 * \brief spatial parameters for the sequential 2p2c test
 */
#ifndef TEST_2P2C_SPATIALPARAMS_HH
#define TEST_2P2C_SPATIALPARAMS_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/2p2c/sequential/properties.hh>
#include <dumux/material/spatialparams/sequentialfv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>

namespace Dumux
{
//forward declaration
template<class TypeTag>
class Test2P2CSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
namespace TTag {
struct Test2P2CSpatialParams {};
}

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Test2P2CSpatialParams> { using type = Test2P2CSpatialParams<TypeTag>; };

}

/*!
 * \ingroup SequentialTwoPTwoCTests
 * \brief spatial parameters for the sequential 2p2c test
 */
template<class TypeTag>
class Test2P2CSpatialParams : public SequentialFVSpatialParams<TypeTag>
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    enum { dim = GridView::dimension };
    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using FieldMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    using PcKrSwCurve = FluidMatrix::LinearMaterialDefault<Scalar>;

public:

    const FieldMatrix& intrinsicPermeability (const Element& element) const
    {
        return constPermeability_;
    }

    double porosity(const Element& element) const
    {
        return 0.2;
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

    Test2P2CSpatialParams(const Problem& problem)
    : SequentialFVSpatialParams<TypeTag>(problem)
    , constPermeability_(0)
    {
        typename PcKrSwCurve::BasicParams params(0/*pcEntry*/, 10000/*pcMax*/);
        pcKrSwCurve_ = std::make_unique<PcKrSwCurve>(params);

        for(int i = 0; i < dim; i++)
        {
            constPermeability_[i][i] = 1e-12;
        }
    }

private:
    FieldMatrix constPermeability_;
    std::unique_ptr<const PcKrSwCurve> pcKrSwCurve_;

};

} // end namespace
#endif
