// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// ## The main program (`main.cc`)
// This file contains the main program flow. In this example, we solve a stationary
// and rotationally symmetric single-phase problem for a sequence of refined grids
// and compute the convergence rates.
// [[content]]
// ### Includes
// [[details]] includes
// [[codeblock]]
#include <config.h>

#include <iostream>
#include <vector>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh> // for GetPropType
#include <dumux/common/parameters.hh> // for getParam
#include <dumux/common/integrate.hh>  // for integrateL2Error

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/pdesolver.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include "properties.hh"
// [[/codeblock]]
// [[/details]]
//
// ### Beginning of the main function
// [[codeblock]]
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // We parse the command line arguments.
    Parameters::init(argc, argv);

    // Convenience alias for the type tag of the problem.
    using TypeTag = Properties::TTag::OnePRotSym;
    // [[/codeblock]]

    // ### Create the grid and the grid geometry
    // [[codeblock]]
    // The grid manager can be used to create a grid from the input file
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();

    // We compute on the leaf grid view.
    const auto& leafGridView = gridManager.grid().leafGridView();

    // instantiate the grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    // [[/codeblock]]

    // ### Initialize the problem and grid variables
    // [[codeblock]]
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // We define a function to update the discrete analytical solution vector
    // using the exactSolution() function in the problem
    const auto updateAnalyticalSolution = [&](auto& pExact, auto& vExact)
    {
        pExact.resize(gridGeometry->numDofs());
        vExact.resize(gridGeometry->elementMapper().size());

        for (const auto& element : elements(gridGeometry->gridView()))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                pExact[scv.dofIndex()] = problem->exactSolution(scv.dofPosition());
            }
            const auto eIdx = gridGeometry->elementMapper().index(element);
            vExact[eIdx] = problem->exactVelocity(element.geometry().center());
        }
    };

    // instantiate and initialize the discrete and exact solution vectors
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector p(gridGeometry->numDofs());
    SolutionVector pExact;
    std::vector<double> vExact;
    updateAnalyticalSolution(pExact, vExact);

    // instantiate and initialize the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(p);
    // [[/codeblock]]

    // ### Initialize VTK output
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, p, problem->name());
    GetPropType<TypeTag, Properties::IOFields>::initOutputModule(vtkWriter);
    vtkWriter.addField(pExact, "pExact"); // add the exact solution to the output fields
    vtkWriter.addField(vExact, "vExact"); // add the exact velocity to the output fields

    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));

    // ### Instantiate the solver
    // We use the `LinearPDESolver` class, which is instantiated on the basis
    // of an assembler and a linear solver. When the `solve` function of the
    // `LinearPDESolver` is called, it uses the assembler and linear
    // solver classes to assemble and solve the linear system around the provided
    // solution and stores the result therein.
    // [[codeblock]]
    using Assembler = FVAssembler<TypeTag, DiffMethod::analytic>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = ILUBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>,
                                               LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
    LinearPDESolver<Assembler, LinearSolver> solver(assembler,  linearSolver);
    solver.setVerbose(false); // suppress output during solve()
    // [[/codeblock]]

    // ### Solution of the problem and error computation
    // The problem is solved by calling `solve` on the instance of `LinearPDESolver`
    // that we have created above. In the following piece of code, we solve the
    // problem on the initial refinement and compute the corresponding L2 error.
    // For a convenient way of computing the L2 error, the function `integrateL2Error`
    // can be used.
    // [[codeblock]]
    solver.solve(p);

    // container to store the L2 errors for the different refinements
    const int numRefinements = getParam<int>("Grid.RefinementSteps");
    std::vector<double> l2Errors(numRefinements);

    // use third order error integration
    constexpr int orderQuadratureRule = 3;

    // compute initial L2 error
    l2Errors[0] = integrateL2Error(*gridGeometry, p, pExact, orderQuadratureRule);
    // [[/codeblock]]

    // This procedure is now repeated for the number of refinements as specified
    // in the input file.
    // [[codeblock]]
    for (int stepIdx = 1; stepIdx < numRefinements; stepIdx++)
    {
        // Globally refine the grid once
        gridManager.grid().globalRefine(1);

        // update the grid geometry, the grid variables and
        // the solution vectors now that the grid has been refined
        gridGeometry->update(gridManager.grid().leafGridView());
        gridVariables->updateAfterGridAdaption(p);

        // this recreates the linear system, i.e. the sizes of
        // the right hand side vector and the Jacobian matrix,
        // and its sparsity pattern.
        assembler->updateAfterGridAdaption();

        p.resize(gridGeometry->numDofs());
        updateAnalyticalSolution(pExact, vExact);

        // solve problem on refined grid
        solver.solve(p);
        // [[/codeblock]]
        // #### Post-processing and output
        // At the end of each refinement step, the convergence
        // rate is printed to the terminal.
        // [[codeblock]]
        // Calculate the L2 error using the numerical solution
        l2Errors[stepIdx] = integrateL2Error(*gridGeometry, p, pExact, orderQuadratureRule);

        // Print the error and convergence rate
        const auto rate = std::log(l2Errors[stepIdx]/l2Errors[stepIdx-1])/std::log(0.5);
        const auto numDofs = gridGeometry->numDofs();
        std::cout << std::setprecision(8) << std::scientific
                  << "-- L2 error for " << std::setw(5) << numDofs << " dofs: " << l2Errors[stepIdx]
                  << ", rate: " << rate
                  << std::endl;
    }
    // [[/codeblock]]

    // After the last refinement, we write the solution to VTK file format on the
    // finest grid and exit the main function.
    // [[codeblock]]
    vtkWriter.write(0.0);

    // program end, return with 0 exit code (success)
    return 0;
}
// [[/codeblock]]
// ### Exception handling
// In this part of the main file we catch and print possible exceptions that could
// occur during the simulation.
// [[details]] error handler
catch (const Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (const Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (const Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
// [[/details]]
// [[/content]]
