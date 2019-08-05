// ## The main file
// This is the main file for the tracer problem. In this file two problems are setup and solved sequentially, first the 1p problem and afterwards the tracer problem. The result of the 1p problem is the pressure distribution in the problem domain. This is used to calculate the volume fluxes, which act as an input for the tracer problem. Based on this volume fluxes the transport of a tracer is calculated in the following tracer problem.

// ### Includes
#include <config.h>

// Both, the problem_1p.hh and the problem_tracer.hh have to be included in the main file.
#include "problem_1p.hh"
#include "problem_tracer.hh"

// Standard header file for C++, to get time and date information.
#include <ctime>

// Standard header file for C++, for in- and output.
#include <iostream>

// Dumux is based on DUNE, the Distributed and Unified Numerics Environment, which provides several grid managers and linear solvers. So we need some includes from that.
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>

// In Dumux a property system is used to specify the model. For this different properties are defined containing type definitions, values and methods. All properties are declared in the file properties.hh.
#include <dumux/common/properties.hh>
// The following file contains the parameter class, which manages the definition of input parameters by a default value, the inputfile or the command line.
#include <dumux/common/parameters.hh>
// The file dumuxmessage.hh contains the class defining the start and end message of the simulation.
#include <dumux/common/dumuxmessage.hh>

// Defines Dumux sequential linear solver backends.
#include <dumux/linear/seqsolverbackend.hh>
// The following file contains the class, which assembles the linear systems for finite volume schemes (box-scheme, tpfa-approximation, mpfa-approximation).
#include <dumux/assembly/fvassembler.hh>
// The containing class in the following file defines the different differentiation methods used to compute the derivatives of the residual.
#include <dumux/assembly/diffmethod.hh>

// The following class is needed to simplify the writing of dumux simulation data to VTK format.
#include <dumux/io/vtkoutputmodule.hh>
// The gridmanager constructs a grid from the information in the input or grid file. There is a specification for the different supported grid managers.
#include <dumux/io/grid/gridmanager.hh>

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // The type tags for the two problems are defined here. They are created in the individual problem files.
    using OnePTypeTag = Properties::TTag::IncompressibleTest;
    using TracerTypeTag = Properties::TTag::TracerTestCC;

    //! initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    //! print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    //! parse command line arguments
    Parameters::init(argc, argv);

    // ### Create the grid

    // A gridmanager tries to create the grid either from a grid file or the input file. Both problems are solved on the same grid. Hence, the grid is only created once using the type tag of the 1p problem.
    GridManager<GetPropType<OnePTypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    //! we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // ### Setup and solving of the 1p problem
    // In the following section the 1p problem is setup and solved. As the result of this problem, the pressure distribution in the problem domain is obtained.

    // #### Setup
    // The finite volume grid geometry, the problem, the linear system, including the jacobian matrix, the residual and the solution vector and the gridvariables are created and initialized.

    // The finite volume geometry builds up the sub control volumes and sub control volume faces for each element of the grid partition.
    using FVGridGeometry = GetPropType<OnePTypeTag, Properties::FVGridGeometry>;
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();

    // In the problem the boundary and initial conditions are defined.
    using OnePProblem = GetPropType<OnePTypeTag, Properties::Problem>;
    auto problemOneP = std::make_shared<OnePProblem>(fvGridGeometry);

    // The jacobian matrix A, the solution vector p and the residual r are parts of the linear system resulting from the Newton's method.
    using JacobianMatrix = GetPropType<OnePTypeTag, Properties::JacobianMatrix>;
    using SolutionVector = GetPropType<OnePTypeTag, Properties::SolutionVector>;
    SolutionVector p(leafGridView.size(0));

    auto A = std::make_shared<JacobianMatrix>();
    auto r = std::make_shared<SolutionVector>();

    // The grid variables store variables on scv and scvf (volume and flux variables).
    using OnePGridVariables = GetPropType<OnePTypeTag, Properties::GridVariables>;
    auto onePGridVariables = std::make_shared<OnePGridVariables>(problemOneP, fvGridGeometry);
    onePGridVariables->init(p);

    // #### Assemling the linear system
    // The assembler is created and inizialized.
    using OnePAssembler = FVAssembler<OnePTypeTag, DiffMethod::analytic>;
    auto assemblerOneP = std::make_shared<OnePAssembler>(problemOneP, fvGridGeometry, onePGridVariables);
    assemblerOneP->setLinearSystem(A, r);

    // The local jacobian and the residual are assembled using Newton's method. The assembly timer stops the time needed for the assembly of the linear system, which is displayed in the terminal output. Further the timer is started to evaluate the total time of the assembly, solving and updating.
    Dune::Timer timer;
    Dune::Timer assemblyTimer; std::cout << "Assembling linear system ..." << std::flush;
    assemblerOneP->assembleJacobianAndResidual(p);
    assemblyTimer.stop(); std::cout << " took " << assemblyTimer.elapsed() << " seconds." << std::endl;

    // We want to solve Ax = -r.
    (*r) *= -1.0;

    // #### Solution
    // For the solution of the linear system the the linear solver "UMFPack" is set as the linear solver. Afterwards the linear system is solved. The time needed to solve the system is recorded by the solverTimer and displayed in the terminal output.
    using LinearSolver = UMFPackBackend;
    Dune::Timer solverTimer; std::cout << "Solving linear system ..." << std::flush;
    auto linearSolver = std::make_shared<LinearSolver>();
    linearSolver->solve(*A, p, *r);
    solverTimer.stop(); std::cout << " took " << solverTimer.elapsed() << " seconds." << std::endl;

    // #### Update and output
    // The grid variables are updated with the new solution.
    Dune::Timer updateTimer; std::cout << "Updating variables ..." << std::flush;
    onePGridVariables->update(p);
    updateTimer.elapsed(); std::cout << " took " << updateTimer.elapsed() << std::endl;


    // We initialize the vtkoutput. Each model has a predefined model specific output with relevant parameters for that model.
    using GridView = GetPropType<OnePTypeTag, Properties::GridView>;
    Dune::VTKWriter<GridView> onepWriter(leafGridView);
    onepWriter.addCellData(p, "p");
    const auto& k = problemOneP->spatialParams().getKField();
    onepWriter.addCellData(k, "permeability");
    onepWriter.write("1p");

    // The timer is stopped and the total time of the simulation as well as the cumulative CPU time is displayed.
    timer.stop();

    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";


    // ### Computation of the volume fluxes
    // The results of the 1p problem are used to calculate the the volume fluxes in the model domain.

    using Scalar =  GetPropType<OnePTypeTag, Properties::Scalar>;
    std::vector<Scalar> volumeFlux(fvGridGeometry->numScvf(), 0.0);

    using FluxVariables =  GetPropType<OnePTypeTag, Properties::FluxVariables>;
    auto upwindTerm = [](const auto& volVars) { return volVars.mobility(0); };

    // We iterate over all elements
    for (const auto& element : elements(leafGridView))
    {
        auto fvGeometry = localView(*fvGridGeometry);
        fvGeometry.bind(element);

        auto elemVolVars = localView(onePGridVariables->curGridVolVars());
        elemVolVars.bind(element, fvGeometry, p);

        auto elemFluxVars = localView(onePGridVariables->gridFluxVarsCache());
        elemFluxVars.bind(element, fvGeometry, elemVolVars);

        // The volume flux is calculated for every subcontrolvolume face, which is not on a Neumann boundary (is not on the boundary or is on a Dirichlet boundary).

        for (const auto& scvf : scvfs(fvGeometry))
        {
            const auto idx = scvf.index();

            if (!scvf.boundary())
            {
                FluxVariables fluxVars;
                fluxVars.init(*problemOneP, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                volumeFlux[idx] = fluxVars.advectiveFlux(0, upwindTerm);
            }
            else
            {
                const auto bcTypes = problemOneP->boundaryTypes(element, scvf);
                if (bcTypes.hasOnlyDirichlet())
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*problemOneP, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                    volumeFlux[idx] = fluxVars.advectiveFlux(0, upwindTerm);
                }
            }
        }
    }


    // ### Setup and solving of the tracer problem

    // #### Setup

    // Similar to the 1p problem first the problem is created and initialized.
    using TracerProblem = GetPropType<TracerTypeTag, Properties::Problem>;
    auto tracerProblem = std::make_shared<TracerProblem>(fvGridGeometry);

    // The volume fluxes calculated in the previous section are used as input for the tracer model.
    tracerProblem->spatialParams().setVolumeFlux(volumeFlux);

    // The solution vector is created and initialized. As the tracer problem is Transient, the initial solution defined in the problem is applied to the solution vector.
    SolutionVector x(leafGridView.size(0));
    tracerProblem->applyInitialSolution(x);
    auto xOld = x;

    //
    using GridVariables = GetPropType<TracerTypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(tracerProblem, fvGridGeometry);
    gridVariables->init(x);

    // Some time loop parameters are read from the input file. The parameter tEnd defines the duration of the simulation, dt the initial time step size and maxDt the maximal time step size. The time step size is adjusted after each time step depending on !!! .
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");

    // The time loop is instantiated.
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // The assembler is created and inizialized with time loop for instationary problem.
    using TracerAssembler = FVAssembler<TracerTypeTag, DiffMethod::analytic, /*implicit=*/false>;
    auto assembler = std::make_shared<TracerAssembler>(tracerProblem, fvGridGeometry, gridVariables, timeLoop);
    assembler->setLinearSystem(A, r);

    // The vtk output module is initialized.
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, tracerProblem->name());
    using IOFields = GetPropType<TracerTypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter);
    using VelocityOutput = GetPropType<TracerTypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    vtkWriter.write(0.0);


    // A check point is set for the time loop.
    timeLoop->setPeriodicCheckPoint(tEnd/10.0);

    // The time loop is started. A new time step is calculated as long as the tEnd is reached. In every single time step the problem is assembled and solved.
    timeLoop->start(); do
    {
        // First the old solution is defined as the solution of the previous time step for storage evaluations.
        assembler->setPreviousSolution(xOld);

        // Then the linear system is assembled.
        Dune::Timer assembleTimer;
        assembler->assembleJacobianAndResidual(x);
        assembleTimer.stop();

        // The linear system A(xOld-xNew) = r is solved.
        Dune::Timer solveTimer;
        SolutionVector xDelta(x);
        linearSolver->solve(*A, xDelta, *r);
        solveTimer.stop();

        // The actual solution is calculated and updated in the grid variables.
        updateTimer.reset();
        x -= xDelta;
        gridVariables->update(x);
        updateTimer.stop();

        // The statistics of the actual time step are displayed.
        const auto elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed();
        std::cout << "Assemble/solve/update time: "
                  <<  assembleTimer.elapsed() << "(" << 100*assembleTimer.elapsed()/elapsedTot << "%)/"
                  <<  solveTimer.elapsed() << "(" << 100*solveTimer.elapsed()/elapsedTot << "%)/"
                  <<  updateTimer.elapsed() << "(" << 100*updateTimer.elapsed()/elapsedTot << "%)"
                  <<  std::endl;

        // The new solution is defined as the old solution.
        xOld = x;
        gridVariables->advanceTimeStep();

        // The time loop is advanced to the next time step.
        timeLoop->advanceTimeStep();

        // Vtk output on check points is written.
        if (timeLoop->isCheckPoint())
            vtkWriter.write(timeLoop->time());

        // Statistics of this time step are reported.
        timeLoop->reportTimeStep();

        // The time step size dt of the next time step is set.
        timeLoop->setTimeStepSize(dt);

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());


    // ### Final Output

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

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
