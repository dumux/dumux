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

#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/integrate.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/io/format.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // We parse the command line arguments.
    Parameters::init(argc, argv);

    // Initialize timer
    Dune::Timer timer;

    // Convenience alias for the type tag of the problem.
    using TypeTag = Properties::TTag::RichardsAnnulus;

    // The grid manager can be used to create a grid from the input file
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();

    // We compute on the leaf grid view.
    const auto& leafGridView = gridManager.grid().leafGridView();

    // instantiate the grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // Initialize the problem and grid variables
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // We define a function to update the discrete analytical solution vector
    // using the exactSolution() function in the problem
    const auto updateAnalyticalSolution = [&](auto& pExact)
    {
        pExact.resize(gridGeometry->numDofs());
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
                pExact[scv.dofIndex()] = problem->exactSolution(scv.dofPosition());
        }
    };

    // instantiate and initialize the discrete and exact solution vectors
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector p(gridGeometry->numDofs()); problem->applyInitialSolution(p);
    SolutionVector pExact; updateAnalyticalSolution(pExact);

    // instantiate and initialize the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(p);

    // initialize VTK output
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, p, problem->name());
    GetPropType<TypeTag, Properties::IOFields>::initOutputModule(vtkWriter);
    vtkWriter.addField(pExact, "pExact"); // add the exact solution to the output fields
    vtkWriter.addVolumeVariable([&](const auto& v){
        return (v.pressure(0) - problem->nonwettingReferencePressure())*100/(v.density(0)*9.81); },
        "head in cm"
    );
    vtkWriter.write(-1.0);

    // instantiate the solver
    using Assembler = FVAssembler<TypeTag, DiffMethod::analytic>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();
    NewtonSolver<Assembler, LinearSolver> solver(assembler, linearSolver);

    // Solution of the problem and error computation
    solver.solve(p);
    vtkWriter.write(0.0);

    const auto initialPressure = problem->headToPressure(getParam<double>("Problem.InitialPressureHeadInCm", -100));
    const auto initialSaturation = problem->saturation(initialPressure);
    std::cout << "Initial uniform profile with saturation: " << initialSaturation << std::endl;
    std::cout << "Time until this profile is reached: "
                  << problem->computeTime(p, *gridVariables, initialSaturation)/(24*60*60)
                  << " days" << std::endl;

    // container to store the L2 errors for the different refinements
    const int numRefinements = getParam<int>("Grid.RefinementSteps");
    std::vector<double> l2Errors(numRefinements);

    // use third order error integration
    constexpr int orderQuadratureRule = 2;

    // compute initial L2 error
    l2Errors[0] = integrateL2Error(*gridGeometry, p, pExact, orderQuadratureRule);
    for (int stepIdx = 1; stepIdx < numRefinements; stepIdx++)
    {
        // Globally refine the grid once
        gridManager.grid().globalRefine(1);

        // update the grid geometry, the grid variables and
        // the solution vectors now that the grid has been refined
        gridGeometry->update(gridManager.grid().leafGridView());
        p.resize(gridGeometry->numDofs());
        problem->applyInitialSolution(p);
        updateAnalyticalSolution(pExact);
        gridVariables->updateAfterGridAdaption(p);

        // this recreates the linear system, i.e. the sizes of
        // the right hand side vector and the Jacobian matrix,
        // and its sparsity pattern.
        assembler->updateAfterGridAdaption();

        // solve problem on refined grid
        solver.solve(p);
        vtkWriter.write(stepIdx);

        std::cout << "Time until this profile is reached: "
                  << problem->computeTime(p, *gridVariables, initialSaturation)/(24*60*60)
                  << " days" << std::endl;

        // Calculate the L2 error using the numerical solution
        l2Errors[stepIdx] = integrateL2Error(*gridGeometry, p, pExact, orderQuadratureRule);

        // Print the error and convergence rate
        const auto rate = std::log(l2Errors[stepIdx]/l2Errors[stepIdx-1])/std::log(0.5);
        const auto numDofs = gridGeometry->numDofs();
        std::cout << Fmt::format(
            "L2 error: {}, rate: {} ({:>5} dofs)\n",
            l2Errors[stepIdx], rate, numDofs
        );

        if (rate < 1.95)
            DUNE_THROW(Dune::Exception, "Convergence rate too small (< 1.95): " << rate);
    }

    std::cout << "Simulation took " << timer.elapsed() << " second." << std::endl;

    return 0;
}
