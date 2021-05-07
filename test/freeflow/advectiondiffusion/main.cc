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
/*!
 * \file
 * \ingroup AdvectionDiffusionTests
 * \brief Test for the advection diffusion model.
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/pdesolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/staggeredfvassembler.hh>

#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    //! define the type tag for this problem
    using FlowTypeTag = Properties::TTag::FLOWTYPETAG;
    using TransportTypeTag = Properties::TTag::TRANSPORTTYPETAG;

    //! initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    //! print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse the command line arguments and input file
    Parameters::init(argc, argv);

    // The hostgrid
    constexpr int dim = 2;
    using Scalar = GetPropType<FlowTypeTag, Properties::Scalar>;
    using GlobalPosition = Dune::FieldVector<Scalar, dim>;
    using HostGrid = Dune::YaspGrid<2, Dune::TensorProductCoordinates<Scalar, dim> >;
    using SubGrid = Dune::SubGrid<dim, HostGrid>;
    using HostGridManager = GridManager<HostGrid>;
    HostGridManager hostGridManager;
    hostGridManager.init();
    auto& hostGrid = hostGridManager.grid();
    Scalar centerCircleRadius = getParam<Scalar>("Problem.CenterCircleRadius");
    GlobalPosition circleCenter = getParam<GlobalPosition>("Problem.CircleCenter");

    // The subgrid selector
    auto gridSelector = [&](const auto& element)
    {
        const auto x = element.geometry().center()[0];
        const auto y = element.geometry().center()[1];
        return std::hypot(x-circleCenter[0], y-circleCenter[1]) > centerCircleRadius;
    };

    Dumux::GridManager<SubGrid> subGridManager;
    subGridManager.init(hostGrid, gridSelector);

    // we compute on the leaf grid view
    const auto& leafGridView = subGridManager.grid().leafGridView();

    ////////////////////////////////////
    // Set and Solve the Flow Problem //
    ////////////////////////////////////

    //! create the finite volume grid geometry
    using FlowGridGeometry = GetPropType<FlowTypeTag, Properties::GridGeometry>;
    auto flowGridGeometry = std::make_shared<FlowGridGeometry>(leafGridView);
    flowGridGeometry->update();

    // the problem (initial and boundary conditions)
    using FlowProblem = GetPropType<FlowTypeTag, Properties::Problem>;
    auto flowProblem = std::make_shared<FlowProblem>(flowGridGeometry);

    // the solution vector
    using FlowSolutionVector = GetPropType<FlowTypeTag, Properties::SolutionVector>;
    FlowSolutionVector xFlow;
    xFlow[FlowGridGeometry::cellCenterIdx()].resize(flowGridGeometry->numCellCenterDofs());
    xFlow[FlowGridGeometry::faceIdx()].resize(flowGridGeometry->numFaceDofs());

    // the grid variables
    using FlowGridVariables = GetPropType<FlowTypeTag, Properties::GridVariables>;
    auto flowGridVariables = std::make_shared<FlowGridVariables>(flowProblem, flowGridGeometry);
    flowGridVariables->init(xFlow);

    // intialize the vtk output module
    StaggeredVtkOutputModule<FlowGridVariables, FlowSolutionVector> flowVtkWriter(*flowGridVariables, xFlow, flowProblem->name() + "_Flow");
    using FlowIOFields = GetPropType<FlowTypeTag, Properties::IOFields>;
    FlowIOFields::initOutputModule(flowVtkWriter); // Add model specific output fields
    flowVtkWriter.write(0.0);

    // the assembler with time loop for instationary problem
    using FlowAssembler = StaggeredFVAssembler<FlowTypeTag, DiffMethod::numeric>;
    auto flowAssembler = std::make_shared<FlowAssembler>(flowProblem, flowGridGeometry, flowGridVariables);

    // the linear solver
    using FlowLinearSolver = Dumux::UMFPackBackend;
    auto flowLinearSolver = std::make_shared<FlowLinearSolver>();

    // the non-linear solver
    using FlowNewtonSolver = Dumux::NewtonSolver<FlowAssembler, FlowLinearSolver>;
    FlowNewtonSolver nonLinearSolver(flowAssembler, flowLinearSolver);

    // linearize & solve
    Dune::Timer timer;
    nonLinearSolver.solve(xFlow);

    // write vtk output
    flowVtkWriter.write(1.0);
    timer.stop();

    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    std::cout << "Flow Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

    /////////////////////////////////////////////////////////////////////////
    /// Pass the Velocity field to the Transport Problem's Spatial Params ///
    /////////////////////////////////////////////////////////////////////////

    using GridIndexType = FlowGridGeometry::GridView::IndexSet::IndexType;
    using SolVectorIndex = FlowSolutionVector::size_type;

    // Iterate through the intersections via the elements to map one dofIdx to each scvf
    GridIndexType intersectionIdx = 0;
    std::map<GridIndexType, SolVectorIndex> scvfToSolVectorIdxMap;
    for (const auto& element : elements(leafGridView))
    {
        for (const auto& intersection : intersections(leafGridView, element))
        {
            // Fill the Map: map key is the SCVF index, map value is the dofIdx
            const SolVectorIndex dofIdx = leafGridView.indexSet().subIndex(intersection.inside(), intersection.indexInInside(), 1);
            scvfToSolVectorIdxMap.insert(std::pair<GridIndexType, SolVectorIndex>(intersectionIdx, dofIdx));
            intersectionIdx++;
        }
    }

    //! create the finite volume grid geometry for the transport problem
    using TransportGridGeometry = GetPropType<TransportTypeTag, Properties::GridGeometry>;
    auto transportGridGeometry = std::make_shared<TransportGridGeometry>(leafGridView);
    transportGridGeometry->update();

    // get Flow problem Solution
    const auto velocitySolution = xFlow[FlowGridGeometry::faceIdx()];
    // Open a vector for the velocity solution at each face, initialize to zero
    std::vector<GlobalPosition> velocityAtScvf(transportGridGeometry->numScvf(), GlobalPosition(0.0));

    // Iterate through each element to get reach each scvf in the transport problem
    for (const auto& element : elements(leafGridView))
    {
        auto fvGeometry = localView(*transportGridGeometry);
        fvGeometry.bind(element);

        // For each scvf we store a velocity vector.
        for (const auto& scvf : scvfs(fvGeometry))
        {
            for (size_t i = 0; i < scvf.unitOuterNormal().size(); i++)
            {
                // The unit Outer Normal is directional, as is the velocity, so we take the absolute value of the unit vector
                using std::abs;
                velocityAtScvf[scvf.index()][i] = velocitySolution[scvfToSolVectorIdxMap[scvf.index()]] * abs(scvf.unitOuterNormal()[i]);
            }
        }
    }

    //! the transport problem (initial and boundary conditions)
    using TransportProblem = GetPropType<TransportTypeTag, Properties::Problem>;
    auto transportProblem = std::make_shared<TransportProblem>(transportGridGeometry);

    // Pass the velocity vector to the spatial params
    transportProblem->spatialParams().setVelocityAtFace(velocityAtScvf);

    ////////////////////////////////////////////////////////////
    // setup instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    //! the solution vector
    using TransportSolutionVector = GetPropType<TransportTypeTag, Properties::SolutionVector>;
    TransportSolutionVector xTransport(transportGridGeometry->numDofs());
    transportProblem->applyInitialSolution(xTransport);
    auto xTransportOld = xTransport;

    //! the grid variables
    using TransportGridVariables = GetPropType<TransportTypeTag, Properties::GridVariables>;
    auto transportGridVariables = std::make_shared<TransportGridVariables>(transportProblem, transportGridGeometry);
    transportGridVariables->init(xTransport);

    //! get some time loop parameters
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");

    //! instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    //! the assembler with time loop for instationary problem
    using TransportAssembler = FVAssembler<TransportTypeTag, DiffMethod::numeric>;
    auto transportAssembler = std::make_shared<TransportAssembler>(transportProblem,
                                                                   transportGridGeometry,
                                                                   transportGridVariables,
                                                                   timeLoop,
                                                                   xTransportOld);

    //! the linear solver
    using TransportLinearSolver = UMFPackBackend;
    auto transportLinearSolver = std::make_shared<TransportLinearSolver>();

    //! pde solver (assemble, solve, update)
    LinearPDESolver transportSolver(transportAssembler, transportLinearSolver);
    transportAssembler->assembleJacobianAndResidual(xTransport);
    transportSolver.reuseMatrix();

    //! intialize the vtk output module
    VtkOutputModule<TransportGridVariables, TransportSolutionVector> transportVtkWriter(*transportGridVariables,
                                                                                        xTransport,
                                                                                        transportProblem->name() + "_Transport" );
    using TransportIOFields = GetPropType<TransportTypeTag, Properties::IOFields>;
    TransportIOFields::initOutputModule(transportVtkWriter); // Add model specific output fields
    transportVtkWriter.write(0.0);

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // run instationary non-linear simulation
    /////////////////////////////////////////////////////////////////////////////////////////////////

    //! start the time loop
    timeLoop->start();
    while (!timeLoop->finished())
    {
        // assemble & solve & update
        transportSolver.solve(xTransport);

        // make the new solution the old solution
        xTransportOld = xTransport;
        transportGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // write vtk output on check points
        transportVtkWriter.write(timeLoop->time());

        // set new dt
        timeLoop->setTimeStepSize(dt);
    }

    timeLoop->finalize(leafGridView.comm());

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    //! print dumux end message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;

}
