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
 *
 * \brief spatial parameters for the test problem for diffusion models.
 */
#ifndef TEST1_FVCA6_SPATIALPARAMETERS_HH
#define TEST1_FVCA6_SPATIALPARAMETERS_HH

#include <dumux/common/properties.hh>
#include <dumux/material/spatialparams/sequentialfv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TestDiffusionSpatialParams3d;

namespace Properties
{
// The spatial parameters TypeTag
namespace TTag {
struct TestDiffusionSpatialParams3d {};
}

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TestDiffusionSpatialParams3d> { using type = TestDiffusionSpatialParams3d<TypeTag>; };

}

/*!
 * \ingroup IMPETtests
 * \brief spatial parameters for the test problem for diffusion models.
 */
template<class TypeTag>
class TestDiffusionSpatialParams3d: public SequentialFVSpatialParams<TypeTag>
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
    using FieldMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    using PcKrSwCurve = FluidMatrix::LinearMaterialDefault<Scalar>;

public:

    const FieldMatrix& intrinsicPermeabilityAtPos (const GlobalPosition& globalPos) const
    {
        return permeability_;
    }

    double porosity(const Element& element) const
    {
        return 1.0;
    }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     *
     * \param globalPos The global coordinates for the given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        return makeFluidMatrixInteraction(pcKrSwCurve_);
    }

    TestDiffusionSpatialParams3d(const Problem& problem)
    : ParentType(problem)
    , pcKrSwCurve_("SpatialParams")
    , permeability_(0)
    {
        // permeability values
        permeability_[0][0] = permeability_[1][1] = permeability_[2][2] = 1.0;
        permeability_[0][1] = permeability_[1][0] = permeability_[1][2] = permeability_[2][1] = 0.5;
    }

private:
    const PcKrSwCurve pcKrSwCurve_;
    FieldMatrix permeability_;
};

} // end namespace
#endif
