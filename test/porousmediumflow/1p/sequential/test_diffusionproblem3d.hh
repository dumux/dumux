// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
 * \brief test problem for diffusion models from the FVCA6 benchmark.
 */
#ifndef DUMUX_TEST_DIFFUSION_3D_PROBLEM_HH
#define DUMUX_TEST_DIFFUSION_3D_PROBLEM_HH

#include <dune/alugrid/grid.hh>

#include <dumux/material/components/constant.hh>

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/lmethod/3dpressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mimetic/pressureproperties.hh>

#include <dumux/porousmediumflow/2p/sequential/diffusion/problem.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/velocity.hh>

#include <dumux/linear/seqsolverbackend.hh>

#include "test_diffusionspatialparams3d.hh"

namespace Dumux
{
/*!
 * \ingroup IMPETtests
 */
template<class TypeTag>
class TestDiffusion3DProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
// Create new type tags
namespace TTag {
struct DiffusionTest { using InheritsFrom = std::tuple<TestDiffusionSpatialParams3d, SequentialTwoP>; };

// set the types for the 2PFA FV method
struct FVTest { using InheritsFrom = std::tuple<DiffusionTest, FVPressureTwoP>; };

// set the types for the MPFA-L FV method
struct FVMPFAL3DTestTypeTag { using InheritsFrom = std::tuple<DiffusionTest, FvMpfaL3dPressureTwoP>; };

// set the types for the mimetic FD method
struct MimeticTest { using InheritsFrom = std::tuple<DiffusionTest, MimeticPressureTwoP>; };

} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DiffusionTest> { using type = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::DiffusionTest> { using type = TestDiffusion3DProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DiffusionTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

template<class TypeTag>
struct LinearSolver<TypeTag, TTag::DiffusionTest> { using type = UMFPackBackend; };

}

/*!
 * \ingroup SequentialProblems
 *
 * \brief test problem for diffusion models from the FVCA6 benchmark.
 */
template<class TypeTag>
class TestDiffusion3DProblem: public DiffusionProblem2P<TypeTag>
{
    using ParentType = DiffusionProblem2P<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Grid = typename GridView::Grid;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pWIdx = Indices::pressureIdx,
        swIdx = Indices::swIdx,
        pressureEqIdx = Indices::pressureEqIdx,
    };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;
    using ScalarSolution = typename SolutionTypes::ScalarSolution;

    TestDiffusion3DProblem(Grid& grid) :
        ParentType(grid), velocity_(*this)
    { }

    //!for this specific problem: initialize the saturation and afterwards the model
    void init()
    {
        this->variables().initialize();
        for (int i = 0; i < this->gridView().size(0); i++)
        {
            this->variables().cellData(i).setSaturation(wPhaseIdx, 1.0);
            this->variables().cellData(i).setSaturation(nPhaseIdx, 0.0);
        }
        this->model().initialize();
    }

    void calculateFVVelocity()
    {
        velocity_.calculateVelocity();
//        velocity_.addOutputVtkFields(this->resultWriter());
    }

    //! \copydoc ParentType::addOutputVtkFields()
    void addOutputVtkFields()
    {
        ScalarSolution *exactPressure = this->resultWriter().allocateManagedBuffer(this->gridView().size(0));

        for(const auto& element : elements(this->gridView()))
        {
            (*exactPressure)[this->elementMapper().index(element)][0] = exact(element.geometry().center());
        }

        this->resultWriter().attachCellData(*exactPressure, "exact pressure");

        return;
    }

    /*!
    * \name Problem parameters
    */
    // \{

    bool shouldWriteRestartFile() const
    { return false; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return 273.15 + 10; // -> 10°C
    }

    // \}


    //! Returns the reference pressure for evaluation of constitutive relations
    Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
    {
        return 1e5; // -> 10°C
    }

    void sourceAtPos(PrimaryVariables &values,const GlobalPosition& globalPos) const
    {
        values = 0;
        using std::sin;
        using std::cos;
        using std::atan;

        double pi = 4.0*atan(1.0);

        values[wPhaseIdx] = -pi*pi*cos(pi*globalPos[0])*cos(pi*(globalPos[1]+1.0/2.0))*sin(pi*(globalPos[2]+1.0/3.0))+
                            3.0*pi*pi*sin(pi*globalPos[0])*sin(pi*(globalPos[1]+1.0/2.0))*sin(pi*(globalPos[2]+1.0/3.0))-
                            pi*pi*sin(pi*globalPos[0])*cos(pi*(globalPos[1]+1.0/2.0))*cos(pi*(globalPos[2]+1.0/3.0));
    }

    /*!
    * \brief Returns the type of boundary condition.
    *
    *
    * BC for saturation equation can be dirichlet (saturation), neumann (flux), or outflow.
    */
    void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
    {
                bcTypes.setAllDirichlet();
    }

    //! set dirichlet condition  (saturation [-])
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
            values[pWIdx] = exact(globalPos);
            values[swIdx] = 1.0;
    }

    //! set neumann condition for phases (flux, [kg/(m^2 s)])
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values = 0;
    }

    Scalar exact (const GlobalPosition& globalPos) const
    {
        using std::sin;
        using std::atan;

        double pi = 4.0*atan(1.0);

        return (1.0+sin(pi*globalPos[0])*sin(pi*(globalPos[1]+1.0/2.0))*sin(pi*(globalPos[2]+1.0/3.0)));
    }

    Dune::FieldVector<Scalar,dim> exactGrad (const GlobalPosition& globalPos) const
    {
        Dune::FieldVector<Scalar,dim> grad(0);
        using std::sin;
        using std::cos;
        using std::atan;

        double pi = 4.0*atan(1.0);
        grad[0] = pi*cos(pi*globalPos[0])*sin(pi*(globalPos[1]+1.0/2.0))*sin(pi*(globalPos[2]+1.0/3.0));
        grad[1] = pi*sin(pi*globalPos[0])*cos(pi*(globalPos[1]+1.0/2.0))*sin(pi*(globalPos[2]+1.0/3.0));
        grad[2] = pi*sin(pi*globalPos[0])*sin(pi*(globalPos[1]+1.0/2.0))*cos(pi*(globalPos[2]+1.0/3.0));

        return grad;
    }

private:
    FVVelocity<TypeTag, GetPropType<TypeTag, Properties::Velocity> > velocity_;
    static constexpr Scalar eps_ = 1e-4;

};
} //end namespace

#endif
