// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <type_traits>
#include <random>

#include <dumux/common/random.hh>

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/discretization/box.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/pdesolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include <examples/diffusion/model.hh>

namespace Dumux {

template<class TypeTag>
class DiffusionTestProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GlobalPosition = typename GridGeometry::LocalView::Element::Geometry::GlobalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;

public:
    DiffusionTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, const Scalar D)
    : ParentType(gridGeometry)
    , diffusionCoefficient_(D)
    {}

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(1.0); }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    Scalar diffusionCoefficient() const
    { return diffusionCoefficient_; }

private:
    Scalar diffusionCoefficient_;
};

} // end namespace Dumux

namespace Dumux::Properties::TTag {

template<class S, int dim>
struct DiffusionTest
{
    using InheritsFrom = std::tuple<DiffusionModel, BoxModel>;

    using Scalar = S;
    using Grid = Dune::YaspGrid<dim>;

    template<class TypeTag>
    using Problem = DiffusionTestProblem<TypeTag>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

} // end namespace Dumux::Properties::TTag

namespace Dumux {

// solve diffusion and return quantities of interest
// here one: total concentration in domain after time tEnd
template<int dim, class Scalar = double>
class DiffusionSolver
{
    using TypeTag = Properties::TTag::DiffusionTest<Scalar, dim>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = UMFPackIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using Solver = Dumux::LinearPDESolver<Assembler, LinearSolver>;
    using TimeLoop = CheckPointTimeLoop<Scalar>;

public:
    DiffusionSolver(const int level)
    {
        Dune::FieldVector<Scalar, Grid::dimension> upperRight(1.0);
        std::array<int, Grid::dimension> cells; cells.fill(2);
        gridManager.init(upperRight, cells);
        for (int i = 0; i < level; ++i)
            gridManager.grid().globalRefine(1);

        gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
    }

    Scalar solve(const Scalar D, const Scalar tEnd, std::size_t numTimeSteps) const
    {
        auto problem = std::make_shared<Problem>(gridGeometry, D);
        SolutionVector sol;
        problem->applyInitialSolution(sol);
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(sol);

        auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, tEnd/numTimeSteps, tEnd);
        auto oldSol = sol;
        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, oldSol);
        auto linearSolver = std::make_shared<LinearSolver>();
        auto solver = std::make_shared<Solver>(assembler, linearSolver);

        assembler->setLinearSystem();
        assembler->assembleJacobian(sol);
        solver->reuseMatrix();

        timeLoop->start(); do
        {
            // assemble & solve
            solver->solve(sol);

            // make the new solution the old solution
            oldSol = sol;
            gridVariables->advanceTimeStep();

            // advance to the time loop to the next step
            timeLoop->advanceTimeStep();

        } while (!timeLoop->finished());

        double amount = 0.0;
        auto fvGeometry = localView(*gridGeometry);
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
                amount += scv.volume()*sol[scv.dofIndex()][0];
        }

        return amount;
    }

private:
    GridManager<Grid> gridManager;
    std::shared_ptr<GridGeometry> gridGeometry;
};

} // end namespace Dumux
