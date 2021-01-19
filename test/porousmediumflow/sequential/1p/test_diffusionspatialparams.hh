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
#ifndef TEST_DIFFUSION_SPATIALPARAMS_HH
#define TEST_DIFFUSION_SPATIALPARAMS_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/sequential/properties.hh>
#include <dumux/material/spatialparams/sequentialfv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TestDiffusionSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
namespace TTag {
struct TestDiffusionSpatialParams {};
}

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TestDiffusionSpatialParams> { using type = TestDiffusionSpatialParams<TypeTag>; };

}

/*!
 * \ingroup IMPETtests
 * \brief spatial parameters for the test problem for diffusion models.
 */
template<class TypeTag>
class TestDiffusionSpatialParams: public SequentialFVSpatialParams<TypeTag>
{
    using ParentType = SequentialFVSpatialParams<TypeTag>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using IndexSet = typename GridView::IndexSet;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CoordScalar = typename Grid::ctype;
    ///@cond false
    using ScalarSolution = typename GetProp<TypeTag, Properties::SolutionTypes>::ScalarSolution;
    ///@endcond
    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld};
    using Element = typename Grid::Traits::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FieldMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    using PcKrSwCurve = FluidMatrix::LinearMaterialDefault<Scalar>;

public:

    const FieldMatrix& intrinsicPermeability (const Element& element) const
    {
        return permeability_[indexSet_.index(element)];
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
        return makeFluidMatrixInteraction(pcKrSwCurve_);
    }

    void initialize(const double delta)
    {
        delta_ = delta;
        permeability_.resize(gridView_.size(0));

        for(const auto& element : elements(gridView_))
        {
            perm(permeability_[indexSet_.index(element)], element.geometry().center());
        }

    }

    template<class Writer>
    void addOutputVtkFields(Writer& writer)
    {
        ScalarSolution *permXX = writer.allocateManagedBuffer(gridView_.size(0));
        ScalarSolution *permXY = writer.allocateManagedBuffer(gridView_.size(0));
        ScalarSolution *permYY = writer.allocateManagedBuffer(gridView_.size(0));

        for(const auto& element : elements(gridView_))
        {
            int eIdxGlobal = indexSet_.index(element);
            (*permXX)[eIdxGlobal][0] = permeability_[eIdxGlobal][0][0];
            (*permXY)[eIdxGlobal][0] = permeability_[eIdxGlobal][0][1];
            (*permYY)[eIdxGlobal][0] = permeability_[eIdxGlobal][1][1];
        }

        writer.attachCellData(*permXX, "permeability-X");
        writer.attachCellData(*permYY, "permeability-Y");
        writer.attachCellData(*permXY, "permeability-Offdiagonal");

        return;
    }

    TestDiffusionSpatialParams(const Problem& problem)
    : ParentType(problem),gridView_(problem.gridView())
    , indexSet_(problem.gridView().indexSet())
    , pcKrSwCurve_("SpatialParams")
    , permeability_(0)
    {}

private:
    void perm (FieldMatrix& perm, const GlobalPosition& globalPos) const
    {
        double rt = globalPos[0]*globalPos[0]+globalPos[1]*globalPos[1];
        perm[0][0] = (delta_*globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1])/rt;
        perm[0][1] = -(1.0 - delta_)*globalPos[0]*globalPos[1]/rt;
        perm[1][0] = perm[0][1];
        perm[1][1] = (globalPos[0]*globalPos[0] + delta_*globalPos[1]*globalPos[1])/rt;
    }

    const GridView gridView_;
    const IndexSet& indexSet_;
    const PcKrSwCurve pcKrSwCurve_;
    std::vector<FieldMatrix> permeability_;
    double delta_;
};

} // end namespace
#endif
