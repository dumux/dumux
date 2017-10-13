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
 * \brief Todo
 */

#ifndef DUMUX_INSTATIONARYNONLINEAR_START_HH
#define DUMUX_INSTATIONARYNONLINEAR_START_HH

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/nonlinear/newtoncontroller.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/io/vtkoutputmodule.hh>

namespace Dumux
{

  /*!
   * \ingroup Simulation
   * \brief Struct that contains the program flow for the solution of an instationary problem.
   *
   * \note Per default we use numerical differentiation for the assembly of the jacobian matrix
   *       and a fully implicit time integration scheme.
   */
template<class TypeTag, DiffMethod diffMethod = DiffMethod::numeric, bool isImplicit = true>
struct InstationaryNonLinearSimulation
{
    static int start(int argc, char** argv)
    {
        // initialize MPI, finalize is done automatically on exit
        const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

        // print dumux start message
        if (mpiHelper.rank() == 0)
            DumuxMessage::print(/*firstCall=*/true);

        // parse command line arguments and input file
        Parameters::init(argc, argv);

        // try to create a grid (from the given grid file or the input file)
        using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
        GridCreator::makeGrid(Parameters::getTree());
        GridCreator::loadBalance();

        ////////////////////////////////////////////////////////////
        // run instationary non-linear problem on this grid
        ////////////////////////////////////////////////////////////

        // we compute on the leaf grid view
        const auto& leafGridView = GridCreator::grid().leafGridView();

        // create the finite volume grid geometry
        using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
        auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
        fvGridGeometry->update();

        // the problem (initial and boundary conditions)
        using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
        auto problem = std::make_shared<Problem>(fvGridGeometry);

        // the solution vector
        using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
        using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
        static constexpr bool isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox);
        static constexpr int dofCodim = isBox ? GridView::dimension : 0;
        SolutionVector x(leafGridView.size(dofCodim));
        problem->applyInitialSolution(x);
        auto xOld = x;

        // the grid variables
        using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
        auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
        gridVariables->init(x, xOld);

        // get some time loop parameters
        using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
        const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
        const auto maxDivisions = getParam<int>("TimeLoop.MaxTimeStepDivisions");
        const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
        auto dt = getParam<Scalar>("TimeLoop.DtInitial");

        // check if we are about to restart a previously interrupted simulation
        Scalar restartTime = 0;
        if (Parameters::getTree().hasKey("Restart") || Parameters::getTree().hasKey("TimeLoop.Restart"))
            restartTime = getParam<Scalar>("TimeLoop.Restart");

        // intialize the vtk output module
        using VtkOutputFields = typename GET_PROP_TYPE(TypeTag, VtkOutputFields);
        VtkOutputModule<TypeTag> vtkWriter(*problem, *fvGridGeometry, *gridVariables, x, problem->name());
        VtkOutputFields::init(vtkWriter); //! Add model specific output fields
        vtkWriter.write(0.0);

        // instantiate time loop
        auto timeLoop = std::make_shared<TimeLoop<Scalar>>(restartTime, dt, tEnd);
        timeLoop->setMaxTimeStepSize(maxDt);

        // the assembler with time loop for instationary problem
        using Assembler = FVAssembler<TypeTag, diffMethod, isImplicit>;
        auto assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables, timeLoop);

        // the linear solver
        using LinearSolver = typename GET_PROP_TYPE(TypeTag, LinearSolver);
        auto linearSolver = std::make_shared<LinearSolver>();

        // the non-linear solver
        using NewtonController = Dumux::NewtonController<TypeTag>;
        using NewtonMethod = Dumux::NewtonMethod<TypeTag, NewtonController, Assembler, LinearSolver>;
        auto newtonController = std::make_shared<NewtonController>(leafGridView.comm(), timeLoop);
        NewtonMethod nonLinearSolver(newtonController, assembler, linearSolver);

        // time loop
        timeLoop->start(); do
        {
            // set previous solution for storage evaluations
            assembler->setPreviousSolution(xOld);

            // try solving the non-linear system
            for (int i = 0; i < maxDivisions; ++i)
            {
                // linearize & solve
                auto converged = nonLinearSolver.solve(x);

                if (converged)
                    break;

                if (!converged && i == maxDivisions-1)
                    DUNE_THROW(Dune::MathError,
                               "Newton solver didn't converge after "
                               << maxDivisions
                               << " time-step divisions. dt="
                               << timeLoop->timeStepSize()
                               << ".\nThe solutions of the current and the previous time steps "
                               << "have been saved to restart files.");
            }

            // make the new solution the old solution
            xOld = x;
            gridVariables->advanceTimeStep();

            // advance to the time loop to the next step
            timeLoop->advanceTimeStep();

            // write vtk output
            vtkWriter.write(timeLoop->time());

            // report statistics of this time step
            timeLoop->reportTimeStep();

            // set new dt as suggested by newton controller
            timeLoop->setTimeStepSize(newtonController->suggestTimeStepSize(timeLoop->timeStepSize()));

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
};

} // end namespace

#endif
