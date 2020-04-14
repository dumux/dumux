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
 * \brief test problem for diffusion models from the FVCA5 benchmark.
 */
#ifndef DUMUX_TEST_2P_PROBLEM_HH
#define DUMUX_TEST_2P_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dumux/material/components/constant.hh>

#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/omethod/2dpressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mimetic/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/problem.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/velocity.hh>

#include "test_diffusionspatialparams.hh"

namespace Dumux {

/*!
 * \ingroup IMPETtests
 */
template<class TypeTag>
class TestDiffusionProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
//// set the types for the 2PFA FV method
// Create new type tags
namespace TTag {

struct FVVelocity2PTestTypeTag { using InheritsFrom = std::tuple<TestDiffusionSpatialParams, FVPressureTwoP>; };

// set the types for the MPFA-O FV method
struct FVMPFAOVelocity2PTestTypeTag { using InheritsFrom = std::tuple<TestDiffusionSpatialParams, FvMpfaO2dPressureTwoP>; };

// set the types for the mimetic FD method
struct MimeticPressure2PTestTypeTag { using InheritsFrom = std::tuple<TestDiffusionSpatialParams, MimeticPressureTwoP>; };
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::FVVelocity2PTestTypeTag> { using type = TestDiffusionProblem<TypeTag>; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::FVVelocity2PTestTypeTag> { using type = Dune::YaspGrid<2>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FVVelocity2PTestTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

//template<class TypeTag>
// struct LinearSolver<TypeTag, TTag::FVMPFAOVelocity2PTestTypeTag> { using type = ILUnBiCGSTABBackend; };
template<class TypeTag>
struct LinearSolver<TypeTag, TTag::FVMPFAOVelocity2PTestTypeTag> { using type = SSORBiCGSTABBackend; };
template<class TypeTag>
struct Problem<TypeTag, TTag::FVMPFAOVelocity2PTestTypeTag> { using type = TestDiffusionProblem<TypeTag>; };
// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::FVMPFAOVelocity2PTestTypeTag> { using type = Dune::YaspGrid<2>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FVMPFAOVelocity2PTestTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};


template<class TypeTag>
struct Problem<TypeTag, TTag::MimeticPressure2PTestTypeTag> { using type = TestDiffusionProblem<TypeTag>; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::MimeticPressure2PTestTypeTag> { using type = Dune::YaspGrid<2>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::MimeticPressure2PTestTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

} // end namespace Properties

/*!
 * \ingroup SequentialProblems
 *
 * \brief test problem for diffusion models from the FVCA5 benchmark.
 *
 * The problem corresponds to Test 2 of the FVCA5 benchmark
 * session, http://www.latp.univ-mrs.fr/fvca5/benchmark/index.html.
 */
template<class TypeTag>
class TestDiffusionProblem: public DiffusionProblem2P<TypeTag>
{
    using ParentType = DiffusionProblem2P<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Grid = typename GridView::Grid;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using WettingPhase = typename GetProp<TypeTag, Properties::FluidSystem>::WettingPhase;

    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        pwIdx = Indices::pwIdx,
        swIdx = Indices::swIdx
    };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using LocalPosition = Dune::FieldVector<Scalar, dim>;

public:
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;
    using ScalarSolution = typename SolutionTypes::ScalarSolution;

    TestDiffusionProblem(Grid& grid) :
        ParentType(grid), velocity_(*this)
    {
        delta_ = getParam<Scalar>("Problem.Delta", 1e-3);
    }

    //!for this specific problem: initialize the saturation and afterwards the model
    void init()
    {
        this->variables().initialize();
        this->spatialParams().initialize(delta_);
        for (int i = 0; i < this->gridView().size(0); i++)
        {
            this->variables().cellData(i).setSaturation(wPhaseIdx, 1.0);
        }
        this->model().initialize();
        velocity_.initialize();
    }

    /*!
    * \name Problem parameters
    */
    // \{

    bool shouldWriteRestartFile() const
    { return false; }


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

        this->spatialParams().addOutputVtkFields(this->resultWriter());

        return;
    }

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

    void source(PrimaryVariables &values,const Element& element) const
    {
        values = 0;

        values[wPhaseIdx] = integratedSource_(element, 4);
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
            values[pwIdx] = exact(globalPos);
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

        Scalar pi = 4.0*atan(1.0);

        return (sin(pi*globalPos[0])*sin(pi*globalPos[1]));
    }

    Dune::FieldVector<Scalar,dim> exactGrad (const GlobalPosition& globalPos) const
        {
        Dune::FieldVector<Scalar,dim> grad(0);
        using std::sin;
        using std::cos;
        using std::atan;

        Scalar pi = 4.0*atan(1.0);
        grad[0] = pi*cos(pi*globalPos[0])*sin(pi*globalPos[1]);
        grad[1] = pi*cos(pi*globalPos[1])*sin(pi*globalPos[0]);

        return grad;
        }

private:

    Scalar integratedSource_(const Element& element, int integrationPoints) const
    {
        Scalar source = 0.;
        LocalPosition localPos(0.0);
        GlobalPosition globalPos(0.0);
        Scalar halfInterval = 1.0/double(integrationPoints)/2.;
        for (int i = 1; i <= integrationPoints; i++)
        {
            for (int j = 1; j <= integrationPoints; j++)
            {
                localPos[0] = double(i)/double(integrationPoints) - halfInterval;
                localPos[1] = double(j)/double(integrationPoints) - halfInterval;
                globalPos = element.geometry().global(localPos);
                source += 1./(integrationPoints*integrationPoints) * evaluateSource_(globalPos);
            }
        }

        return source;
    }

    Scalar evaluateSource_(const GlobalPosition& globalPos) const
    {
        Scalar temp = temperatureAtPos(globalPos);
        Scalar referencePress = referencePressureAtPos(globalPos);

        using std::sin;
        using std::cos;
        using std::atan;

        Scalar pi = 4.0 * atan(1.0);
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        Scalar dpdx = pi * cos(pi * x) * sin(pi * y);
        Scalar dpdy = pi * sin(pi * x) * cos(pi * y);
        Scalar dppdxx = -pi * pi * sin(pi * x) * sin(pi * y);
        Scalar dppdxy = pi * pi * cos(pi * x) * cos(pi * y);
        Scalar dppdyx = dppdxy;
        Scalar dppdyy = dppdxx;
        Scalar kxx = (delta_* x*x + y*y)/(x*x + y*y);
        Scalar kxy = -(1.0 - delta_) * x * y / (x*x + y*y);
        Scalar kyy = (x*x + delta_*y*y)/(x*x + y*y);
        Scalar dkxxdx = 2 * x * y*y * (delta_ - 1.0)/((x*x + y*y) * (x*x + y*y));
        Scalar dkyydy = 2 * x*x * y * (delta_ - 1.0)/((x*x + y*y) * (x*x + y*y));
        Scalar dkxydx = (1.0 - delta_) * y * (x*x - y*y) /((x*x + y*y) * (x*x + y*y));
        Scalar dkxydy = (1.0 - delta_) * x * (y*y - x*x) /((x*x + y*y) * (x*x + y*y));

        Scalar fx = dkxxdx * dpdx + kxx * dppdxx + dkxydx * dpdy + kxy * dppdyx;
        Scalar fy = dkxydy * dpdx + kxy * dppdxy + dkyydy * dpdy + kyy * dppdyy;

        return -(fx + fy) / WettingPhase::viscosity(temp, referencePress) * WettingPhase::density(temp, referencePress);
    }

    Scalar delta_;
    FVVelocity<TypeTag, GetPropType<TypeTag, Properties::Velocity> > velocity_;
};
} //end namespace

#endif
