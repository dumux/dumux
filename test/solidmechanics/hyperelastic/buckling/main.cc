// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>
#include <iostream>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#if USE_GRIDFORMAT
#include <dumux/io/gridwriter.hh>
#endif

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::HyperelasticityTest;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // initialize parameter tree
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    x = 0.0;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // initialize the vtk output module
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    vtkWriter.addField(x, "d");

    // the assembler with time loop for a transient problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    // the linear solver
    using LAT = LinearAlgebraTraitsFromAssembler<Assembler>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LAT>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver);

    nonLinearSolver->solve(x);
    vtkWriter.write(0.0);

    const auto incrementFraction = 0.01;
    const auto length = gridGeometry->bBoxMax()[0] - gridGeometry->bBoxMin()[0];
    auto increment = incrementFraction * length;
    const auto maxDisplacement = 0.2 * length;

    auto direction = 1.0;

    auto currentLength = length;
    double counter = 0.0;
    auto oldSol = x;
    while (currentLength > length - maxDisplacement)
    {
        bool converged = false;
        int maxReduce = 10;
        while (!converged && maxReduce > 0)
        {
            currentLength = currentLength + direction * increment;
            std::cout << "Current displacement: " << currentLength - length << std::endl;
            problem->setCurrentDisplacement(currentLength-length);

            try {
                nonLinearSolver->solve(x);
                converged = true;
            } catch (Dumux::NumericalProblem& e)
            {
                std::cout << "Retry" << std::endl;
                currentLength = currentLength - direction * increment;
                increment = 0.5*increment;
                --maxReduce;
                x = oldSol;
                gridVariables->update(x);
            }

            if (maxReduce == 0)
                DUNE_THROW(Dumux::NumericalProblem, "Newton didn't converge.");
        }

        if (maxReduce > 9 && increment < incrementFraction * length)
            increment = 2.0*increment;

        counter += 1.0;
        vtkWriter.write(counter);
        oldSol = x;

        if (currentLength >= length + maxDisplacement || currentLength <= length - maxDisplacement)
            direction *= -1.0;
    }

    return 0;
}
