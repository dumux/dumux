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
 * \ingroup TracerTests
 * \brief Test for the 2p tracer CC model
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>

#include <test/porousmediumflow/2p/implicit/incompressible/problem.hh>
#include "problem_tracer.hh"

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/amgbackend.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/method.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

int main(int argc, char** argv)
{
    using namespace Dumux;

    //! define the type tags for this problem
    using TwoPTypeTag = Properties::TTag::TwoPIncompressibleTpfa;
    using TracerTypeTag = Properties::TTag::TwoPTracerTestTpfa;
    //! initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    //! print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TwoPTypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // get some time loop parameters
    using Scalar =  GetPropType<TwoPTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = getParam<Scalar>("Restart.Time", 0);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    ////////////////////////////////////////////////////////////
    // set 2p Problem
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TwoPTypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // the problem (initial and boundary conditions)
    using TwoPProblem = GetPropType<TwoPTypeTag, Properties::Problem>;
    auto twoPProblem = std::make_shared<TwoPProblem>(gridGeometry);

    // the solution vector
    using TwoPSolutionVector = GetPropType<TwoPTypeTag, Properties::SolutionVector>;
    TwoPSolutionVector p(gridGeometry->numDofs());
    twoPProblem->applyInitialSolution(p);
    auto pOld = p;

    // maybe update the interface parameters
    if constexpr (ENABLEINTERFACESOLVER)
        twoPProblem->spatialParams().updateMaterialInterfaces(p);

    // the grid variables
    using TwoPGridVariables = GetPropType<TwoPTypeTag, Properties::GridVariables>;
    auto twoPGridVariables = std::make_shared<TwoPGridVariables>(twoPProblem, gridGeometry);
    twoPGridVariables->init(p);

    // intialize the vtk output module
    using TwoPIOFields = GetPropType<TwoPTypeTag, Properties::IOFields>;

    // use non-conforming output for the test with interface solver
    const auto ncOutput = getParam<bool>("Problem.UseNonConformingOutput", false);
    VtkOutputModule<TwoPGridVariables, TwoPSolutionVector> twoPVtkWriter(*twoPGridVariables, p, twoPProblem->name()+"_2p", "",
                                                                        (ncOutput ? Dune::VTK::nonconforming : Dune::VTK::conforming));

    TwoPIOFields::initOutputModule(twoPVtkWriter); //!< Add model specific output fields
    twoPVtkWriter.write(0.0);

    // the assembler with time loop for instationary problem
    using TwoPAssembler = FVAssembler<TwoPTypeTag, DiffMethod::numeric>;
    auto twoPAssembler = std::make_shared<TwoPAssembler>(twoPProblem, gridGeometry, twoPGridVariables, timeLoop, pOld);

    // the linear solver
    using TwoPLinearSolver = AMGBiCGSTABBackend<LinearSolverTraits<GridGeometry>>;
    auto twoPLinearSolver = std::make_shared<TwoPLinearSolver>(leafGridView, gridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<TwoPAssembler, TwoPLinearSolver>;
    NewtonSolver nonLinearSolver(twoPAssembler, twoPLinearSolver);

    ////////////////////////////////////////////////////////////
    // set tracer Problem
    ////////////////////////////////////////////////////////////

    //! the problem (initial and boundary conditions)
    using TracerProblem = GetPropType<TracerTypeTag, Properties::Problem>;
    auto tracerProblem = std::make_shared<TracerProblem>(gridGeometry);

    //! the solution vector
    using TracerSolutionVector = GetPropType<TracerTypeTag, Properties::SolutionVector>;
    TracerSolutionVector x(leafGridView.size(0));
    tracerProblem->applyInitialSolution(x);
    auto xOld = x;

    //! initialize the flux, density and saturation vectors
    std::vector<Scalar> volumeFlux_(gridGeometry->numScvf(), 0.0);
    std::vector<Scalar> density_(gridGeometry->numScv(), 0.0);
    std::vector<Scalar> saturation_(gridGeometry->numScv(), 0.0);

    //! the grid variables
    using TracerGridVariables = GetPropType<TracerTypeTag, Properties::GridVariables>;
    auto tracerGridVariables = std::make_shared<TracerGridVariables>(tracerProblem, gridGeometry);
    tracerGridVariables->init(x);

    // the linear solver
    using TracerLinearSolver = AMGBiCGSTABBackend<LinearSolverTraits<GridGeometry>>;
    auto tracerLinearSolver = std::make_shared<TracerLinearSolver>(leafGridView, gridGeometry->dofMapper());

     //! the linear system
    using JacobianMatrix = GetPropType<TracerTypeTag, Properties::JacobianMatrix>;
    auto A = std::make_shared<JacobianMatrix>();
    auto r = std::make_shared<TracerSolutionVector>();

    //! the assembler with time loop for instationary problem
    using TracerAssembler = FVAssembler<TracerTypeTag, DiffMethod::analytic, /*implicit=*/false>;
    auto tracerAssembler = std::make_shared<TracerAssembler>(tracerProblem, gridGeometry, tracerGridVariables, timeLoop, xOld);
    tracerAssembler->setLinearSystem(A, r);

    // set the flux, density and saturation from the 2p problem
    tracerProblem->spatialParams().setVolumeFlux(volumeFlux_);
    tracerProblem->spatialParams().setDensity(density_);
    tracerProblem->spatialParams().setSaturation(saturation_);

    //! initialize the vtk output module
    VtkOutputModule<TracerGridVariables, TracerSolutionVector> vtkWriter(*tracerGridVariables, x, tracerProblem->name());
    using TracerIOFields = GetPropType<TracerTypeTag, Properties::IOFields>;
    TracerIOFields::initOutputModule(vtkWriter); //!< Add model specific output fields
    using VelocityOutput = GetPropType<TracerTypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*tracerGridVariables));
    vtkWriter.write(0.0);

    //! set some check points for the time loop
    timeLoop->setPeriodicCheckPoint(tEnd/50.0);

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem
    ////////////////////////////////////////////////////////////

    // time loop
    timeLoop->start(); do
    {
        // solve the non-linear system with time step control
        nonLinearSolver.solve(p, *timeLoop);

        // make the new solution the old solution
        pOld = p;
        twoPGridVariables->advanceTimeStep();

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by the Newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));


        // loop over elements to compute fluxes, saturations, densities for tracer
        using FluxVariables = GetPropType<TwoPTypeTag, Properties::FluxVariables>;
        auto upwindTerm = [](const auto& volVars) { return volVars.mobility(0); };
        for (const auto& element : elements(leafGridView))
        {
            ////////////////////////////////////////////////////////////
            // compute volume fluxes for the tracer model
            ///////////////////////////////////////////////////////////

            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bind(element);

            auto elemVolVars = localView(twoPGridVariables->curGridVolVars());
            elemVolVars.bind(element, fvGeometry, p);

            auto elemFluxVars = localView(twoPGridVariables->gridFluxVarsCache());
            elemFluxVars.bind(element, fvGeometry, elemVolVars);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                const auto idx = scvf.index();

                if (!scvf.boundary())
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*twoPProblem, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                            volumeFlux_[idx] = fluxVars.advectiveFlux(0, upwindTerm);
                }
                else
                {
                    const auto bcTypes = twoPProblem->boundaryTypes(element, scvf);
                    if (bcTypes.hasOnlyDirichlet())
                    {
                        FluxVariables fluxVars;
                        fluxVars.init(*twoPProblem, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                                volumeFlux_[idx] = fluxVars.advectiveFlux(0, upwindTerm);
                    }
                }
            }

            ////////////////////////////////////////////////////////////
            // compute densities and saturations for the tracer model
            ///////////////////////////////////////////////////////////
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                const auto idx = scv.dofIndex();

                density_[idx] = volVars.density(0);
                saturation_[idx] = volVars.saturation(0);
            }
        }

        ////////////////////////////////////////////////////////////
        // solve tracer problem on the same grid
        ////////////////////////////////////////////////////////////

        // set the flux from the 2p problem
        tracerProblem->spatialParams().setVolumeFlux(volumeFlux_);
        tracerProblem->spatialParams().setDensity(density_);
        tracerProblem->spatialParams().setSaturation(saturation_);

        Dune::Timer tracerAssembleTimer;
        tracerAssembler->assembleJacobianAndResidual(x);
        tracerAssembleTimer.stop();

        // solve the linear system A(xOld-xNew) = r
        Dune::Timer solveTimer;
        TracerSolutionVector xDelta(x);
        tracerLinearSolver->solve(*A, xDelta, *r);
        solveTimer.stop();

        // update solution and grid variables
        Dune::Timer updateTimer;
        updateTimer.reset();
        x -= xDelta;
        tracerGridVariables->update(x);
        updateTimer.stop();

        // statistics
        Dune::Timer assembleTimer;
        const auto elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed();
        std::cout << "Assemble/solve/update time: "
                  <<  assembleTimer.elapsed() << "(" << 100*assembleTimer.elapsed()/elapsedTot << "%)/"
                  <<  solveTimer.elapsed() << "(" << 100*solveTimer.elapsed()/elapsedTot << "%)/"
                  <<  updateTimer.elapsed() << "(" << 100*updateTimer.elapsed()/elapsedTot << "%)"
                  <<  std::endl;

        // make the new solution the old solution
        xOld = x;
        tracerGridVariables->advanceTimeStep();

        // advance the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        twoPVtkWriter.write(timeLoop->time());
        vtkWriter.write(timeLoop->time());

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    //! print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
