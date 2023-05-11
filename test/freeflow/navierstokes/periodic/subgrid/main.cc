// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Periodic Pore Scale Simulations
 * \ingroup Basic Pore Structure
 * \brief A test for the evaluation of a pore geometry
 */

#include <config.h>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/gmshwriter.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_sp.hh>
#include <dumux/io/grid/gridmanager_sub.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/vtk/function.hh>
#include <dumux/io/vtk/intersectionwriter.hh>
#include <dumux/freeflow/navierstokes/velocityoutput.hh>

#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectionentityset.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/linear/pdesolver.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>

#include "properties.hh"

#include <dumux/common/metadata.hh>

template<class GlobalPosition>
std::string positionToString(const GlobalPosition& pos, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    for (const auto& coord : pos)
        out << std::fixed << coord;
    return out.str();
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using FlowMomentumTypeTag = Properties::TTag::FlowMomentumModel;
    using FlowMassTypeTag = Properties::TTag::FlowMassModel;

    // Initialize MPI, finalize is done automatically on exit. Print dumux start message
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    ///////////////////////////////////////////
    // Step 1: Build the Pore Scale Geometry //
    ///////////////////////////////////////////

    // Define dimension and type containers
    constexpr int dim = 2;
    using Scalar = GetPropType<FlowMomentumTypeTag, Properties::Scalar>;
    using IntArray = std::array<int, dim>;
    using GlobalPosition = Dune::FieldVector<Scalar, dim>;

    // Collect the image type and number
    const std::string imgFileName = getParam<std::string>("Grid.PBMFileName");
    const IntArray hostCells = getParam<IntArray>("Grid.Cells");

    ////// Set up the host-grid
    using HostGrid = Dune::SPGrid<Scalar, dim>;
    using HostGridManager = GridManager<HostGrid>;
    HostGridManager hostGridManager;
    hostGridManager.init("");
    auto& hostGrid = hostGridManager.grid();
    const auto& hostGridView = hostGrid.leafGridView();
    using SingleHostGridGeometry = CCTpfaFVGridGeometry<typename HostGrid::LeafGridView>;
    auto singleHostGridGeometry = std::make_shared<SingleHostGridGeometry>(hostGridView);

    // Set up the Sub-grid
    std::cout << "Constructing SubGrid from binary image ("
              << imgFileName << ") and a SP hostgrid" << std::endl;

    using SubGrid = Dune::SubGrid<dim, HostGrid>;
    using SubGridManager = Dumux::GridManager<SubGrid>;
    SubGridManager subGridManager;
    const auto img = NetPBMReader::readPBM(imgFileName);
    const bool marked = getParam<bool>("Grid.Marker", 0);

    auto gridSelector = [&hostGrid, &img, &marked](const auto& element)
    { return img[hostGrid.leafGridView().indexSet().index(element)] == marked; };

    subGridManager.init(hostGrid, gridSelector);
    const auto& poreScaleGridView = subGridManager.grid().leafGridView();

    //////////////////////////////////////////////
    // Step 2: Calculate Pore Scale Flow Field  //
    //////////////////////////////////////////////

    // run the stationary non-linear problem on the pore scale grid
    std::cout << "\nSolve the flow problem on the pore scale geometry \n";

    // create the finite volume grid geometry
    using PoreScaleMomentumGridGeometry = GetPropType<FlowMomentumTypeTag, Properties::GridGeometry>;
    auto poreScaleMomentumGridGeometry = std::make_shared<PoreScaleMomentumGridGeometry>(poreScaleGridView);
    using PoreScaleMassGridGeometry = GetPropType<FlowMassTypeTag, Properties::GridGeometry>;
    auto poreScaleMassGridGeometry = std::make_shared<PoreScaleMassGridGeometry>(poreScaleGridView);
    std::cout << "periodicity subgrid mass: " << std::boolalpha << poreScaleMassGridGeometry->isPeriodic() << "\n";
    std::cout << "periodicity subgrid momentum: " << std::boolalpha << poreScaleMomentumGridGeometry->isPeriodic() << "\n";

    // Figure out the full porosity
    double fullPorosity = double(poreScaleMassGridGeometry->numScv()) / double(singleHostGridGeometry->numScv());
    std::cout << "Porosity of the subgrid: " << fullPorosity << "\n";

    // the coupling manager
    using CouplingManager = GetPropType<FlowMomentumTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problem (boundary conditions)
    using MomentumProblem = GetPropType<FlowMomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(poreScaleMomentumGridGeometry, couplingManager);
    using MassProblem = GetPropType<FlowMassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(poreScaleMassGridGeometry, couplingManager);

    // the solution vector
    using Traits = MultiDomainTraits<FlowMomentumTypeTag, FlowMassTypeTag>;
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector flowSolutionVector;
    flowSolutionVector[momentumIdx].resize(poreScaleMomentumGridGeometry->numDofs());
    flowSolutionVector[massIdx].resize(poreScaleMassGridGeometry->numDofs());

    // the grid variables
    using MomentumGridVariables = GetPropType<FlowMomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, poreScaleMomentumGridGeometry);
    using MassGridVariables = GetPropType<FlowMassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, poreScaleMassGridGeometry);

    // initialize coupling stencil first and then grid variables (need coupling variables)
    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), flowSolutionVector);
    massGridVariables->init(flowSolutionVector[massIdx]);
    momentumGridVariables->init(flowSolutionVector[momentumIdx]);

    // initialize the vtk output module
    using IOFields = GetPropType<FlowMassTypeTag, Properties::IOFields>;
    VtkOutputModule flowVtkWriter(*massGridVariables, flowSolutionVector[massIdx], massProblem->name());
    IOFields::initOutputModule(flowVtkWriter); // Add model specific output fields
    flowVtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());

    // the assembler with time loop for instationary problem
    using FlowAssembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto flowAssembler = std::make_shared<FlowAssembler>(std::make_tuple(momentumProblem, massProblem),
                                                         std::make_tuple(poreScaleMomentumGridGeometry, poreScaleMassGridGeometry),
                                                         std::make_tuple(momentumGridVariables, massGridVariables),
                                                         couplingManager);

    // the linear solver
    using FlowLinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<FlowAssembler>>;
    auto flowLinearSolver = std::make_shared<FlowLinearSolver>();

    // the non-linear solver
    using FlowNewtonSolver = MultiDomainNewtonSolver<FlowAssembler, FlowLinearSolver, CouplingManager>;
    FlowNewtonSolver flowNonLinearSolver(flowAssembler, flowLinearSolver, couplingManager);

    // linearize & solve
    Dune::Timer flowtimer;
    flowNonLinearSolver.solve(flowSolutionVector);

    // write vtk output
    flowVtkWriter.write(0.0);
    flowtimer.stop();

    const auto& flowComm = Dune::MPIHelper::getCommunication();
    std::cout << "Solve took " << flowtimer.elapsed() << " seconds on "
            << flowComm.size() << " processes.\n"
            << "The cumulative CPU time for the flow problem was " << flowtimer.elapsed()*flowComm.size() << " seconds.\n";
    std::cout << "Flow problem on the pore scale geometry is solved. \n";

    //////////////////////////////////////////////
    // Step 3: Extend to a tiled grid geometry  //
    //////////////////////////////////////////////

    // Set up the host-grid
    HostGridManager tiledHostGridManager;
    IntArray tiledCells = hostCells;
    IntArray numRepeats = getParam<IntArray>("Grid.NumRepeats");
    GlobalPosition lowerLeft = getParam<GlobalPosition>("Grid.LowerLeft");
    GlobalPosition upperRight = getParam<GlobalPosition>("Grid.UpperRight");
    for (int i = 0; i < dim; i++)
    {
        tiledCells[i] *= numRepeats[i];
        lowerLeft[i] -= (numRepeats[i]-1)/2;
        upperRight[i] += (numRepeats[i]-1)/2;
    }
    std::bitset<dim> periodic = getParam<std::bitset<dim>>("Grid.Periodic");
    tiledHostGridManager.init(lowerLeft, upperRight, tiledCells, "", 1, periodic);
    auto& tiledHostGrid = tiledHostGridManager.grid();

    SubGridManager tiledSubGridManager;

    // rewrite img file according to stamp request
    std::vector<bool> stampedImg;
    const int numRows = img.header().nRows;
    const int numCols = img.header().nCols;

    // Reorder the single img file to fill the stamped vector
    for (int stampY = 0; stampY < numRepeats[1]; stampY++)
        for (int j = 0; j < numRows; j++)
            for (int stampX = 0; stampX < numRepeats[0]; stampX++)
                for (int i = 0; i < numCols; i++)
                    stampedImg.push_back(img[ (j * numCols + i) ]);

    // Create the element selector for a repeated image
    auto repeatedElementSelector = [&tiledHostGrid, &stampedImg, &marked](const auto& element)
    { return stampedImg[tiledHostGrid.leafGridView().indexSet().index(element) ] == marked; };

    tiledSubGridManager.init(tiledHostGrid, repeatedElementSelector, "");
    const auto& tiledPoreScaleGridView = tiledSubGridManager.grid().leafGridView();
    using TiledPoreGridGeometry = GetPropType<FlowMassTypeTag, Properties::GridGeometry>;
    auto tiledPoreGridGeometry = std::make_shared<TiledPoreGridGeometry>(tiledPoreScaleGridView);
    using TiledPoreMomentumGridGeometry = GetPropType<FlowMomentumTypeTag, Properties::GridGeometry>;
    auto tiledPoreMomentumGridGeometry = std::make_shared<TiledPoreMomentumGridGeometry>(tiledPoreScaleGridView);
    std::cout << "periodicity tiledPoreGridGeometry: " << std::boolalpha << tiledPoreGridGeometry->isPeriodic() << "\n";

    //! initialize the vtk output module for the Scales
    VtkOutputModuleBase<TiledPoreGridGeometry> tiledVtkWriter(*tiledPoreGridGeometry, "test_ff_periodic_subgrid_tiled", "Double");

    // Extend the face velocity solution to the tiled grid
    using VelocityVector = GlobalPosition;
    std::vector<VelocityVector> faceVelocityVectorTiled(poreScaleMomentumGridGeometry->numScv(), VelocityVector(0.0));
    std::vector<GlobalPosition> positionSet(poreScaleMomentumGridGeometry->numScv(), GlobalPosition(0.0));

    for (const auto& element : elements(poreScaleMomentumGridGeometry->gridView()))
    {
        auto fvGeometry = localView(*poreScaleMomentumGridGeometry);
        fvGeometry.bind(element);
        auto elemVolVars = localView(momentumGridVariables->curGridVolVars());
        elemVolVars.bind(element, fvGeometry, flowSolutionVector);

        for (const auto& scv : scvs(fvGeometry))
        {
            positionSet[scv.index()] = scv.dofPosition();
            faceVelocityVectorTiled[scv.index()][scv.dofAxis()] = elemVolVars[scv].velocity();
        }
    }

    // Create an unordered map to store the solution values associated with x and y positions
    std::unordered_map<std::string, VelocityVector> velocitySolutionMap;
    for (int i = 0; i < positionSet.size(); i++)
    {
        std::string positionString = positionToString(positionSet[i]);
        velocitySolutionMap.insert(std::make_pair(positionString, faceVelocityVectorTiled[i]));
    }

    // Use the map to get the velocity onto the tiled faces, map solutions from the unit cell.
    GlobalPosition unitCellBBoxMax = singleHostGridGeometry->bBoxMax();
    GlobalPosition unitCellBBoxMin = singleHostGridGeometry->bBoxMin();
    std::vector<VelocityVector> tiledVelocityFace(tiledPoreGridGeometry->numScvf(), VelocityVector(-50.0));
    for (const auto& element : elements(tiledPoreGridGeometry->gridView()))
    {
        auto fvGeometry = localView(*tiledPoreGridGeometry);
        fvGeometry.bind(element);
        for (const auto& scvf : scvfs(fvGeometry))
        {
            // reduce the position back to the original
            GlobalPosition globalPos = scvf.center();
            GlobalPosition outerNormal = scvf.unitOuterNormal();
            for (int axis = 0; axis < dim; axis++)
            {
                while ( globalPos[axis] > (unitCellBBoxMax[axis] + (1e-10 * outerNormal[axis])) )
                    globalPos[axis] -= 1.0;
            }

            for (int axis = 0; axis < dim; axis++)
            {
                while ( globalPos[axis] < (unitCellBBoxMin[axis] + (1e-10 * outerNormal[axis])) )
                    globalPos[axis] += 1.0;
            }

            std::string positionString = positionToString(globalPos);
            auto faceSol = velocitySolutionMap.find(positionString);
            if (faceSol != velocitySolutionMap.end())
                tiledVelocityFace[scvf.index()] = faceSol->second;
            else
            {
                std::cout << "No solution value found for (" << globalPos[0] << ", " << globalPos[1] << ")"
                          << "with unitNormal Vector " << scvf.unitOuterNormal()[0] << ", " << scvf.unitOuterNormal()[1] << ") \n";
            }
        }
    }

    // calculate cell-center-averaged velocities
    std::vector<VelocityVector> tiledVelocityVectorCC(tiledPoreGridGeometry->elementMapper().size(), VelocityVector(0.0));
    for (const auto& element : elements(tiledPoreGridGeometry->gridView()))
    {
        auto fvGeometry = localView(*tiledPoreGridGeometry);
        fvGeometry.bindElement(element);
        unsigned int elementIdx = tiledPoreGridGeometry->elementMapper().index(element);

        // calculate velocities
        VelocityVector velocityTemp(0.0);
        for (auto&& scvf : scvfs(fvGeometry))
        {
            const int dofIdxFace = scvf.index();
            const auto numericalSolutionFace = tiledVelocityFace[dofIdxFace];
            velocityTemp += numericalSolutionFace;
        }

        for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            tiledVelocityVectorCC[elementIdx][dimIdx] = velocityTemp[dimIdx] * 0.5; // faces are equidistant to cell center
    }

    tiledVtkWriter.addField(tiledVelocityVectorCC, "velocity", Vtk::Precision::float64);
    tiledVtkWriter.write(1.0);

    Parameters::print();

    return 0;

} // end main
