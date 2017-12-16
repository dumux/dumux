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
 * \brief spatial parameters for the test problem for diffusion models.
 */
#ifndef TEST_DIFFUSION_SPATIALPARAMS_HH
#define TEST_DIFFUSION_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/sequential/properties.hh>
#include <dumux/material/spatialparams/sequentialfv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TestDiffusionSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(TestDiffusionSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TestDiffusionSpatialParams, SpatialParams, TestDiffusionSpatialParams<TypeTag>);

// Set the material law
SET_PROP(TestDiffusionSpatialParams, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef LinearMaterial<Scalar> RawMaterialLaw;
public:
    typedef EffToAbsLaw<RawMaterialLaw> type;
};
}

/*!
 * \ingroup IMPETtests
 * \brief spatial parameters for the test problem for diffusion models.
 */
template<class TypeTag>
class TestDiffusionSpatialParams: public SequentialFVSpatialParams<TypeTag>
{
    typedef SequentialFVSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    ///@cond false
    typedef typename GET_PROP(TypeTag, SolutionTypes)::ScalarSolution ScalarSolution;
    ///@endcond
    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    const FieldMatrix& intrinsicPermeability (const Element& element) const
    {
        return permeability_[indexSet_.index(element)];
    }

    double porosity(const Element& element) const
    {
        return 0.2;
    }


    // return the parameter object for the Brooks-Corey material law which depends on the position
    const MaterialLawParams& materialLawParams(const Element &element) const
    {
            return materialLawParams_;
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

    TestDiffusionSpatialParams(const GridView& gridView)
    : ParentType(gridView),gridView_(gridView), indexSet_(gridView.indexSet()), permeability_(0)
    {
        // residual saturations
        materialLawParams_.setSwr(0.0);
        materialLawParams_.setSnr(0.0);

        // parameters for the linear entry pressure function
        materialLawParams_.setEntryPc(0);
        materialLawParams_.setMaxPc(0);
    }

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
    MaterialLawParams materialLawParams_;
    std::vector<FieldMatrix> permeability_;
    double delta_;
};

} // end namespace
#endif
