// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <config.h>
#include <array>
#include <iostream>
#include <type_traits>

#include <dune/common/exceptions.hh>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_ug.hh>

#include <test/freeflow/shallowwater/errors.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // the property tag
    using TypeTag = Properties::TTag::TYPETAG;

    // create grid
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // problem, in which we define the boundary and initial conditions.
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    problem->applyInitialSolution(x);
    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // We then initialize the predefined model-specific output vtk output.
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables,x, problem->name());
    vtkWriter.addField(problem->getExactWaterDepth(), "exactWaterDepth");
    vtkWriter.addField(problem->getExactVelocityX(), "exactVelocityX");
    vtkWriter.addField(problem->getExactVelocityY(), "exactVelocityY");
    problem->updateAnalyticalSolution();
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(0.0);

    // assembler & solver
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = LINEARSOLVER<LinearSolverTraits<GridGeometry>,
                                      LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = [&]{
        if constexpr (std::is_constructible_v<LinearSolver, decltype(leafGridView), decltype(gridGeometry->dofMapper())>)
            return std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());
        else
            return std::make_shared<LinearSolver>();
    }();

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    nonLinearSolver.solve(x);
    vtkWriter.write(1.0);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto errors = ShallowWaterTest::computeErrors(*gridGeometry, x, [&](const auto& element){
        const auto eIdx = gridGeometry->elementMapper().index(element);
        return std::array<Scalar, 3>{{
            problem->getExactWaterDepth()[eIdx], problem->getExactVelocityX()[eIdx], problem->getExactVelocityY()[eIdx]
        }};
    });
    if (mpiHelper.rank() == 0)
        ShallowWaterTest::printErrors(errors);

    const auto errorTolerance = getParam<Scalar>("Problem.ErrorTolerance", -1.0);
    if (errorTolerance > 0.0 && !(errors.l1[0] < errorTolerance && errors.l1[1] < errorTolerance))
        DUNE_THROW(Dune::Exception, "L1 errors of water depth and velocity exceed " << errorTolerance);

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
