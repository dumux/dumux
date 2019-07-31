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
 * \brief test for the 1p2c box model
 */
#include <config.h>
#include <ctime>
#include <iostream>

#include <dune/common/float_cmp.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include "1p2csaltinstabilityproblem.hh"
#include <dumux/discretization/method.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

#include <dumux/discretization/evalgradients.hh>

void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
                    errorMessageOut += progName;
                    errorMessageOut += " [options]\n";
                    errorMessageOut += errorMsg;
                    errorMessageOut += "\n\nThe List of Mandatory arguments for this program is:\n"
                                        "\t-tEnd                               The end of the simulation [s] \n"
                                        "\t-dtInitial                          The initial timestep size [s] \n"
                                        "\t-gridFile                           The file name of the file containing the grid \n"
                                        "\t                                        definition in DGF format\n"
                                        "\t-SpatialParameters.lensLowerLeftX   Dimension of the lens [m] \n"
                                        "\t-SpatialParameters.lensLowerLeftY   Dimension of the lens [m] \n"
                                        "\t-SpatialParameters.lensUpperRightX  Dimension of the lens [m] \n"
                                        "\t-SpatialParameters.lensUpperRighty  Dimension of the lens [m] \n"
                                        "\n";

        std::cout << errorMessageOut
                  << "\n";
    }
}

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::SaltInstability1p2cProblemTypeTag;

    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // initialize parameter tree
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(fvGridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(fvGridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");

    // intialize the vtk output module
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter); //!< Add model specific output fields
    vtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables, timeLoop);

    // the linear solver
    using LinearSolver = ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);

    // GridView, getting dimWorld, Indices
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    static constexpr int dimWorld = GridView::dimensionworld;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    // time loop
    timeLoop->start(); do
    {
       // begin: required for dispersion example
        const auto numElements = fvGridGeometry->elementMapper().size();
        using FV = Dune::FieldMatrix<Scalar,dimWorld,dimWorld>;
        std::vector<FV> dispTensor(numElements, FV(0.0));

        // first time step without dispersion
        if (timeLoop->time()>0)
        {
        // velocity field and Disp.coeff calculation
        for (const auto& element : elements(leafGridView))
        {
          auto fvGeometry = localView(*fvGridGeometry);
          fvGeometry.bind(element);

          auto elemVolVars = localView(gridVariables->prevGridVolVars());
          elemVolVars.bind(element, fvGeometry, x);

          auto elemSol = elementSolution(element,elemVolVars,fvGeometry);

          auto gradP = evalGradients(element,element.geometry(),fvGeometry.fvGridGeometry(),elemSol,element.geometry().center())[Indices::pressureIdx];

          //interpolate the solution
          const auto& localBasis = fvGridGeometry->feCache().get(element.geometry().type()).localBasis();

          // evaluate the shape functions at the scv center
          const auto localPos = element.geometry().local(element.geometry().center());
          std::vector< Dune::FieldVector<Scalar, 1> > shapeValues;
          localBasis.evaluateFunction(localPos, shapeValues);

           Scalar rho(0.0);
           for (auto&& scv : scvs(fvGeometry))
           {
              const auto& volVars = elemVolVars[scv];
              rho += volVars.density(0)*shapeValues[scv.indexInElement()][0];
           }

           gradP.axpy(-rho, problem->gravity());

           auto velocity = gradP;
           velocity *= problem->spatialParams().permeabilityAtPos(localPos);

           Scalar mu(0.0);
           for (auto&& scv : scvs(fvGeometry))
           {
              const auto& volVars = elemVolVars[scv];
              mu += volVars.viscosity(0)*shapeValues[scv.indexInElement()][0];
           }

           velocity /= mu;

           /////////////////////////////////
           //calculate dispersion tensor
           static const auto alphaL = getParam<Scalar>("Problem.AlphaL");
           static const auto alphaT = getParam<Scalar>("Problem.AlphaT");

           std::array<Scalar,2> dispersivity;
           dispersivity[0] = alphaL;
           dispersivity[1] = alphaT;

           const auto eIdx = fvGridGeometry->elementMapper().index(element);

           //matrix multiplication of the velocity at the interface: vv^T
           Scalar velocityProduct;
           for (int i=0; i < dimWorld; i++)
               for (int j = 0; j < dimWorld; j++){
                   velocityProduct = velocity[i]*velocity[j];
                   // if (velocityProduct < 0) velocityProduct *= -1;
                   dispTensor[eIdx][i][j] = velocityProduct;}

           //normalize velocity product --> vv^T/||v||, [m/s]
           Scalar vNorm = velocity.two_norm();

           dispTensor[eIdx] /= vNorm;
           if (vNorm < 1e-20)
               dispTensor[eIdx] = 0;

           //multiply with dispersivity difference: vv^T/||v||*(alphaL - alphaT), [m^2/s] --> alphaL = longitudinal disp., alphaT = transverse disp.
           dispTensor[eIdx] *= (dispersivity[0] - dispersivity[1]);

           //add ||v||*alphaT to the main diagonal:vv^T/||v||*(alphaL - alphaT) + ||v||*alphaT, [m^2/s]
           for (int i = 0; i < dimWorld; i++)
               dispTensor[eIdx][i][i] += vNorm*dispersivity[1];

        } // end dispersionTensor calculation
        problem->spatialParams().setDispTensors(dispTensor);

        // set in spatial params
        }

        // end: required for dispersion example


        // set previous solution for storage evaluations
        assembler->setPreviousSolution(xOld);

        // solve the non-linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by the newton solver
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

catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
              "). Most likely, the DGF file name is wrong "
              "or the DGF file is corrupted, "
              "e.g. missing hash at end of file or wrong number (dimensions) of entries."
              << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
