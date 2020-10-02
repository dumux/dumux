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
 * \ingroup ShallowWaterTests
 * \brief A test for the shallow water model (bowl).
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/io/grid/gridmanager.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/amgbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>

#include "properties.hh"

template<class SolutionVector, int i>
class SolutionComponent
{
    const SolutionVector& sol_;
public:
    SolutionComponent(const SolutionVector& sol) : sol_(sol) {}
    auto operator[] (std::size_t idx) const { return sol_[idx][i]; }
    auto size() const { return sol_.size(); }
};

//! Compute the analytical solution at time t
template<class SolutionVector, class Problem>
SolutionVector computeAnalyticalSolution(const double t, const Problem& problem)
{
    const auto& gg = problem.gridGeometry();
    SolutionVector exactSolution(gg.numDofs());
    for (const auto& element : elements(gg.gridView()))
    {
        const auto eIdx = gg.elementMapper().index(element);
        const auto& globalPos = element.geometry().center();
        exactSolution[eIdx] = problem.analyticalSolution(t, globalPos);
    }
    return exactSolution;
}

//! Compute L2 error wrt the analytical solution at time t
template<class SolutionVector, class GridGeometry>
typename SolutionVector::block_type
computeL2Error(const double t,
               const SolutionVector& exactSolution,
               const SolutionVector& curSol,
               const GridGeometry& gridGeometry)
{
    typename SolutionVector::block_type l2Error(0.0);
    for (const auto& element : elements(gridGeometry.gridView(), Dune::Partitions::interior))
    {
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bindElement(element);

        using std::pow;
        for (auto&& scv : scvs(fvGeometry))
        {
            auto localDiffSq = exactSolution[scv.dofIndex()] - curSol[scv.dofIndex()];
            for (int i = 0; i < localDiffSq.size(); ++i)
                localDiffSq[i] *= localDiffSq[i];
            l2Error += localDiffSq;
        }
    }
    // sum over processes if we are running with multiple processes in parallel
    if (gridGeometry.gridView().comm().size() > 0)
        l2Error = gridGeometry.gridView().comm().sum(l2Error);

    // square root of sum of squared errors is our absolute discrete l2 error
    for (int i = 0; i < l2Error.size(); ++i)
        l2Error[i] = std::sqrt(l2Error[i]);
    return l2Error;
}

//! Compute and print L2 error wrt the analytical solution at time t
template<class SolutionVector, class GridGeometry>
void computeAndPrintL2Error(const double t,
                            const SolutionVector& exactSolution,
                            const SolutionVector& curSol,
                            const GridGeometry& gridGeometry)
{
    const auto l2Error = computeL2Error(t, exactSolution, curSol, gridGeometry);
    auto numElements = gridGeometry.gridView().size(0);
    // sum over processes if we are running with multiple processes in parallel
    if (gridGeometry.gridView().comm().size() > 0)
        numElements = gridGeometry.gridView().comm().sum(numElements);

    if (gridGeometry.gridView().comm().rank() == 0)
        std::cout << "L2 error in (h, u, v) at t = " <<  t << " seconds for " << numElements << " elements: "
                  << std::scientific << l2Error << std::endl;
}

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::Bowl;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = 2.0*problem->oscillationPeriodInSeconds();
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // intialize the vtk output module and analytical solution
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables,x, problem->name());
    auto exactSolution = computeAnalyticalSolution<SolutionVector>(0.0, *problem);
    const auto exactWaterDepth = SolutionComponent<SolutionVector, 0>(exactSolution);
    const auto exactVelocityX = SolutionComponent<SolutionVector, 1>(exactSolution);
    const auto exactVelocityY = SolutionComponent<SolutionVector, 2>(exactSolution);
    vtkWriter.addField(exactWaterDepth, "exactWaterDepth");
    vtkWriter.addField(exactVelocityX, "exactVelocityX");
    vtkWriter.addField(exactVelocityY, "exactVelocityY");
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = AMGBiCGSTABBackend<LinearSolverTraits<GridGeometry>>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    //! Compute L2 error for the initial solution (should be zero because initial solution is exact)
    computeAndPrintL2Error(0.0, exactSolution, x, *gridGeometry);

    // time loop
    timeLoop->start(); do
    {
        nonLinearSolver.solve(x,*timeLoop);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // update the analytical solution and print current l2 error
        exactSolution = computeAnalyticalSolution<SolutionVector>(timeLoop->time(), *problem);
        computeAndPrintL2Error(timeLoop->time(), exactSolution, x, *gridGeometry);

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output every 10 time steps
        if (!(timeLoop->timeStepIndex() % 10))
            vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton controller
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));


    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

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
