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
 * \ingroup SequentialTwoPTests
 * \brief spatial parameters for the explicit transport test
 */
#ifndef TEST_TRANSPORT_SPATIALPARAMS_HH
#define TEST_TRANSPORT_SPATIALPARAMS_HH

#include <dumux/common/properties.hh>
#include <dumux/material/spatialparams/sequentialfv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>

namespace Dumux
{
/*!
 * \ingroup SequentialTwoPTests
 * \brief spatial parameters for the explicit transport test
 */
//forward declaration
template<class TypeTag>
class TestTransportSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
namespace TTag {
struct TestTransportSpatialParams {};
}

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TestTransportSpatialParams> { using type = TestTransportSpatialParams<TypeTag>; };

}

/*!
 * \ingroup IMPETtests
 * \brief spatial parameters for the explicit transport test
 */
template<class TypeTag>
class TestTransportSpatialParams: public SequentialFVSpatialParams<TypeTag>
{
    using ParentType = SequentialFVSpatialParams<TypeTag>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Element = typename Grid::Traits::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::LinearMaterialDefault<Scalar>;

public:


    Scalar intrinsicPermeability (const Element& element) const
    {
        return 1e-5;
    }

    double porosity(const Element& element) const
    {
        return 0.2;
    }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     * \param globalPos A global coordinate vector
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition &globalPos) const
    {
        return makeFluidMatrixInteraction(pcKrSwCurve_);
    }

    TestTransportSpatialParams(const Problem& problem)
    : ParentType(problem), pcKrSwCurve_("SpatialParams")
    {}

private:
    const PcKrSwCurve pcKrSwCurve_;
};

} // end namespace
#endif
