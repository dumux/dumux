// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief test for the pore network model
 */
#include <config.h>

#include <iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/initialize.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/io/grid/gridmanager_sub.hh>

#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/io/grid/porenetwork/subgriddata.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/porenetwork/common/pnmvtkoutputmodule.hh>
#include <dumux/porenetwork/common/boundaryflux.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////

    // parse command line arguments
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    using HostGridManager = PoreNetwork::GridManager<3>;
    HostGridManager hostGridManager;
    hostGridManager.init();
    const auto hostGridView = hostGridManager.grid().leafGridView();
    const auto hostGridData = hostGridManager.getGridData();

    using SubGrid = Dune::SubGrid<1, typename HostGridManager::Grid>;
    using SubGridManager = GridManager<SubGrid>;

    auto elementSelectorVoid = [&](const auto& element)
    {
        return hostGridData->getParameter(element, "ThroatDomainType") == 0;
    };

    auto elementSelectorSolid = [&](const auto& element)
    {
        return hostGridData->getParameter(element, "ThroatDomainType") == 1;
    };

    // for debugging only
    auto elementSelectorCoupling = [&](const auto& element)
    {
        return hostGridData->getParameter(element, "ThroatDomainType") == 2;
    };

    auto& hostGrid = hostGridManager.grid();
    SubGridManager subgridManagerVoid; subgridManagerVoid.init(hostGrid, elementSelectorVoid);
    SubGridManager subgridManagerSolid; subgridManagerSolid.init(hostGridManager.grid(), elementSelectorSolid);
    SubGridManager subgridManagerCoupling; subgridManagerCoupling.init(hostGridManager.grid(), elementSelectorCoupling);
    const auto voidGridView = subgridManagerVoid.grid().leafGridView();
    const auto solidGridView = subgridManagerSolid.grid().leafGridView();
    const auto couplingGridView = subgridManagerCoupling.grid().leafGridView(); // for debugging only

    std::cout << "Void grid has " << voidGridView.size(0) << " elements" << std::endl;
    std::cout << "Solid grid has " << solidGridView.size(0) << " elements" << std::endl;

    using SubGridData = PoreNetwork::SubGridData<HostGridManager::Grid, SubGridManager::Grid>;
    auto voidGridData = std::make_shared<SubGridData>(subgridManagerVoid.grid(), hostGridData);
    auto solidGridData = std::make_shared<SubGridData>(subgridManagerSolid.grid(), hostGridData);

    using SolidTypeTag = Properties::TTag::PNMSolidModel;
    using VoidTypeTag = Properties::TTag::PNMVoidModel;

    // create the finite volume grid geometry
    using VoidGridGeometry = GetPropType<VoidTypeTag, Properties::GridGeometry>;
    auto voidGridGeometry = std::make_shared<VoidGridGeometry>(voidGridView, *voidGridData);
    voidGridGeometry->update(voidGridView,*voidGridData);

    using SolidGridGeometry = GetPropType<SolidTypeTag, Properties::GridGeometry>;
    auto solidGridGeometry = std::make_shared<SolidGridGeometry>(solidGridView,*solidGridData);
    solidGridGeometry->update(solidGridView,*solidGridData);

    // the spatial parameters
    using VoidSpatialParams = GetPropType<VoidTypeTag, Properties::SpatialParams>;
    auto voidSpatialParams = std::make_shared<VoidSpatialParams>(voidGridGeometry, *voidGridData);

    using SolidSpatialParams = GetPropType<SolidTypeTag, Properties::SpatialParams>;
    auto solidSpatialParams = std::make_shared<SolidSpatialParams>(solidGridGeometry, *solidGridData);

    // the coupling manager
    using Traits = MultiDomainTraits<SolidTypeTag, VoidTypeTag>;
    using CouplingManager = PoreNetwork::PNMHeatTransferCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>();
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector sol;

    // the problem (boundary conditions)
    using VoidProblem = GetPropType<VoidTypeTag, Properties::Problem>;
    auto voidProblem = std::make_shared<VoidProblem>(voidGridGeometry, voidSpatialParams, couplingManager);
    sol[CouplingManager::voidDomainIdx].resize(voidProblem->gridGeometry().numDofs());
    voidProblem->applyInitialSolution(sol[CouplingManager::voidDomainIdx]);

    using SolidProblem = GetPropType<SolidTypeTag, Properties::Problem>;
    auto solidProblem = std::make_shared<SolidProblem>(solidGridGeometry, solidSpatialParams, couplingManager);
    sol[CouplingManager::solidDomainIdx].resize(solidProblem->gridGeometry().numDofs());
    solidProblem->applyInitialSolution(sol[CouplingManager::solidDomainIdx]);

    auto solOld = sol;

    // initialize the coupling manager
    couplingManager->init(solidProblem, voidProblem,
                          hostGridView, *hostGridData,
                          voidGridView, solidGridView, sol);

    // the grid variables
    using VoidGridVariables = GetPropType<VoidTypeTag, Properties::GridVariables>;
    auto voidGridVariables = std::make_shared<VoidGridVariables>(voidProblem, voidGridGeometry);
    voidGridVariables->init(sol[CouplingManager::voidDomainIdx]);

    using SolidGridVariables = GetPropType<SolidTypeTag, Properties::GridVariables>;
    auto solidGridVariables = std::make_shared<SolidGridVariables>(solidProblem, solidGridGeometry);
    solidGridVariables->init(sol[CouplingManager::solidDomainIdx]);

    // initialize the vtk output modules
    using Scalar = typename Traits::Scalar;
    using VoidVtkWriter = PoreNetwork::VtkOutputModule<
        VoidGridVariables, GetPropType<VoidTypeTag, Properties::FluxVariables>, decltype(sol[CouplingManager::voidDomainIdx])
    >;
    VoidVtkWriter voidVtkWriter(*voidGridVariables, sol[CouplingManager::voidDomainIdx],  voidProblem->name());
    GetPropType<VoidTypeTag, Properties::IOFields>::initOutputModule(voidVtkWriter);
    using SolidVtkWriter = PoreNetwork::VtkOutputModule<
        SolidGridVariables, GetPropType<SolidTypeTag, Properties::FluxVariables>, decltype(sol[CouplingManager::solidDomainIdx])
    >;
    SolidVtkWriter solidVtkWriter(*solidGridVariables, sol[CouplingManager::solidDomainIdx],  solidProblem->name());
    GetPropType<SolidTypeTag, Properties::IOFields>::initOutputModule(solidVtkWriter);

    // get some time loop parameters
    using Scalar = typename Traits::Scalar;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd", 1.0);
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize", 1.0);
    auto dt = getParam<Scalar>("TimeLoop.DtInitial", 1.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    static const bool isStationary = getParam<bool>("Problem.IsStationary", true);

    // the assembler
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = isStationary ? std::make_shared<Assembler>(std::make_tuple(solidProblem, voidProblem),
                                                                std::make_tuple(solidGridGeometry,
                                                                                voidGridGeometry),
                                                                std::make_tuple(solidGridVariables,
                                                                                voidGridVariables),
                                                                couplingManager)
                                  : std::make_shared<Assembler>(std::make_tuple(solidProblem, voidProblem),
                                                                std::make_tuple(solidGridGeometry,
                                                                                voidGridGeometry),
                                                                std::make_tuple(solidGridVariables,
                                                                                voidGridVariables),
                                                                couplingManager,
                                                                timeLoop, solOld);

    // the linear solver
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits,
                                               LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // pass gridvars to void problem
    voidProblem->setGridVariables(voidGridVariables);

    // pass gridVars to coupling manager
    couplingManager->setGridVariables(std::make_tuple(solidGridVariables, voidGridVariables));

    // add dof indices and corresponding host grid vertex indices
    std::vector<std::size_t> voidDofIdx(voidGridView.size(1));
    std::vector<std::size_t> voidHostGridVertexIdx(voidGridView.size(1));
    std::vector<Scalar> voidConductionSource(voidGridView.size(1));
    std::vector<Scalar> voidConvectionSource(voidGridView.size(1));
    for (const auto& v : vertices(voidGridView))
    {
        const auto vIdx = voidGridView.indexSet().index(v);
        voidDofIdx[vIdx] = vIdx;
        voidHostGridVertexIdx[vIdx] = hostGridView.indexSet().index(v.impl().hostEntity());
    }

    auto voidFVGeometry = localView(*voidGridGeometry);
    auto voidElemVolVars = localView(voidGridVariables->curGridVolVars());
    for (const auto& element : elements(voidGridView))
    {
        voidFVGeometry.bind(element);
        voidElemVolVars.bind(element, voidFVGeometry, sol[CouplingManager::voidDomainIdx]);
        couplingManager->bindCouplingContext(CouplingManager::voidDomainIdx, element);

        for (const auto& scv : scvs(voidFVGeometry))
        {
            if (couplingManager->isCoupledPore(CouplingManager::voidDomainIdx, scv.dofIndex()))
            {
                voidConductionSource[scv.dofIndex()] += couplingManager->conductionSource(CouplingManager::voidDomainIdx, element,
                                                                                          voidFVGeometry, voidElemVolVars, scv)
                                                                                          * scv.volume();

                voidConvectionSource[scv.dofIndex()] += couplingManager->convectionSource(CouplingManager::voidDomainIdx, element,
                                                                                          voidFVGeometry, voidElemVolVars, scv)
                                                                                          * scv.volume();
            }
        }
    }

    voidVtkWriter.addField(voidDofIdx, "dofIdx", Vtk::FieldType::vertex);
    voidVtkWriter.addField(voidHostGridVertexIdx, "hostGridVertexIdx", Vtk::FieldType::vertex);
    voidVtkWriter.addField(voidConductionSource, "conductionSource", Vtk::FieldType::vertex);
    voidVtkWriter.addField(voidConvectionSource, "convectionSource", Vtk::FieldType::vertex);
    voidVtkWriter.addField(voidGridGeometry->poreVolume(), "poreVolume", Vtk::FieldType::vertex);
    voidVtkWriter.addField(voidGridGeometry->throatCrossSectionalArea(), "throatArea", Vtk::FieldType::element);

    std::vector<std::size_t> solidDofIdx(solidGridView.size(1));
    std::vector<std::size_t> solidHostGridVertexIdx(solidGridView.size(1));
    std::vector<Scalar> solidConductionSource(solidGridView.size(1));
    std::vector<Scalar> solidConvectionSource(solidGridView.size(1));
    for (const auto& v : vertices(solidGridView))
    {
        const auto vIdx = solidGridView.indexSet().index(v);
        solidDofIdx[vIdx] = vIdx;
        solidHostGridVertexIdx[vIdx] = hostGridView.indexSet().index(v.impl().hostEntity());
    }

    auto solidFVGeometry = localView(*solidGridGeometry);
    auto solidElemVolVars = localView(solidGridVariables->curGridVolVars());
    for (const auto& element : elements(solidGridView))
    {
        solidFVGeometry.bind(element);
        solidElemVolVars.bind(element, solidFVGeometry, sol[CouplingManager::solidDomainIdx]);
        couplingManager->bindCouplingContext(CouplingManager::solidDomainIdx, element);

        for (const auto& scv : scvs(solidFVGeometry))
        {
            if (couplingManager->isCoupledPore(CouplingManager::solidDomainIdx, scv.dofIndex()))
            {
                solidConductionSource[scv.dofIndex()] += couplingManager->conductionSource(CouplingManager::solidDomainIdx, element,
                                                                                           solidFVGeometry, solidElemVolVars, scv)
                                                                                           * scv.volume();

                solidConvectionSource[scv.dofIndex()] += couplingManager->convectionSource(CouplingManager::solidDomainIdx, element,
                                                                                           solidFVGeometry, solidElemVolVars, scv)
                                                                                           * scv.volume();
            }
        }
    }

    solidVtkWriter.addField(solidDofIdx, "dofIdx", Vtk::FieldType::vertex);
    solidVtkWriter.addField(solidHostGridVertexIdx, "hostGridVertexIdx", Vtk::FieldType::vertex);
    solidVtkWriter.addField(solidConductionSource, "conductionSource", Vtk::FieldType::vertex);
    solidVtkWriter.addField(solidConvectionSource, "convectionSource", Vtk::FieldType::vertex);

    voidVtkWriter.write(0);
    solidVtkWriter.write(0);

    // for debugging
    Dune::VTKWriter<std::decay_t<decltype(couplingGridView)>> couplingWriter(couplingGridView);
    std::vector<std::size_t> couplingDebugVertexIndices(couplingGridView.size(1));
    std::vector<std::size_t> couplingDebugConnectionIds(couplingGridView.size(0));
    std::vector<std::size_t> connectionIdToDebugElemIdx(couplingGridView.size(0));
    std::vector<Scalar> connectionArea(couplingGridView.size(0));
    std::vector<Scalar> connectionLength(couplingGridView.size(0));
    for (const auto& v : vertices(couplingGridView))
        couplingDebugVertexIndices[couplingGridView.indexSet().index(v)] = hostGridView.indexSet().index(v.impl().hostEntity());

    for (const auto& e : elements(couplingGridView))
    {
        const auto hostGridElementIdx = hostGridView.indexSet().index(e.impl().hostEntity());
        const auto id = couplingManager->couplingMapper().hostGridElementIndexToGlobalId().at(hostGridElementIdx);
        const auto debugElemIdx = couplingGridView.indexSet().index(e);
        couplingDebugConnectionIds[debugElemIdx] = id;
        connectionIdToDebugElemIdx[id] = debugElemIdx;

        connectionArea[debugElemIdx] = couplingManager->couplingMapper().connectionInfo()[id].connectionArea;
        connectionLength[debugElemIdx] = couplingManager->couplingMapper().connectionInfo()[id].connectionLength;
    }
    couplingWriter.addVertexData(couplingDebugVertexIndices, "hostGridVertexIdx");
    couplingWriter.addCellData(connectionArea, "connectionArea");
    couplingWriter.addCellData(connectionLength, "connectionLength");
    couplingWriter.addCellData(couplingDebugConnectionIds, "connectionId");


    std::fill(voidConductionSource.begin(), voidConductionSource.end(), 0.0);
    std::fill(voidConvectionSource.begin(), voidConvectionSource.end(), 0.0);
    std::fill(solidConductionSource.begin(), solidConductionSource.end(), 0.0);
    std::fill(solidConvectionSource.begin(), solidConvectionSource.end(), 0.0);

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    const auto voidBoundaryFlux = PoreNetwork::BoundaryFlux(
        *voidGridVariables, assembler->localResidual(CouplingManager::voidDomainIdx), sol[CouplingManager::voidDomainIdx]
    );
    const auto solidBoundaryFlux = PoreNetwork::BoundaryFlux(
        *solidGridVariables, assembler->localResidual(CouplingManager::solidDomainIdx), sol[CouplingManager::solidDomainIdx]
    );

    if (isStationary)
    {
        // solve the non-linear system without time step control
        nonLinearSolver.solve(sol);

        if (couplingGridView.size(0) != couplingManager->couplingMapper().connectionInfo().size())
            DUNE_THROW(Dune::InvalidStateException, "Wrong number of connections");

        std::vector<Scalar> solidConnectionFlux(couplingGridView.size(0), 0.0);
        std::vector<Scalar> voidConnectionFlux(couplingGridView.size(0), 0.0);
        std::vector<Scalar> tVoidMinusTSolid(couplingGridView.size(0), 0.0);
        std::vector<bool> solidConnectionVisited(couplingGridView.size(0), false);
        std::vector<bool> voidConnectionVisited(couplingGridView.size(0), false);

        Scalar totalSolidInletHeatFlux = 0.0;

        couplingWriter.addCellData(voidConnectionFlux, "connectionFluxVoid");
        couplingWriter.addCellData(solidConnectionFlux, "connectionFluxSolid");
        couplingWriter.addCellData(tVoidMinusTSolid, "tVoidMinusTSolid");

        for (const auto& element : elements(solidGridView))
        {
            solidFVGeometry.bind(element);
            solidElemVolVars.bind(element, solidFVGeometry, sol[CouplingManager::solidDomainIdx]);
            couplingManager->bindCouplingContext(CouplingManager::solidDomainIdx, element);

            for (const auto& scv : scvs(solidFVGeometry))
            {
                totalSolidInletHeatFlux += solidProblem->robinInletHeatFlux(element, solidFVGeometry, solidElemVolVars, scv) * scv.volume();
                if (couplingManager->isCoupledPore(CouplingManager::solidDomainIdx, scv.dofIndex()))
                {
                    // iterate over all connection throats
                    const auto& connections = couplingManager->couplingMapper().solidToVoidConnections(scv.dofIndex());
                    const auto& contextForPore = couplingManager->couplingContext()[scv.dofIndex()];

                    for (const auto& [localConnectionIdx, connection] : enumerate(connections))
                    {
                        if (!solidConnectionVisited[connection.id])
                        {
                            solidConnectionVisited[connection.id] = true;
                            const auto& context = contextForPore[localConnectionIdx];
                            const Scalar tVoid = sol[CouplingManager::voidDomainIdx][connection.voidVertexIdx][1];
                            const Scalar tSolid = sol[CouplingManager::solidDomainIdx][connection.solidVertexIdx][0];

                            tVoidMinusTSolid[connection.id] = tVoid - tSolid;

                            if (voidProblem->sourceMode() == decltype(voidProblem->sourceMode())::convection)
                                solidConnectionFlux[connection.id] = couplingManager->convectiveHeatFluxForOneConnection(connection, context);
                            else
                            {
                                const Scalar t = couplingManager->getConnectionTransmissiblity(CouplingManager::solidDomainIdx, connection, solidElemVolVars, scv);
                                const Scalar deltaT = tSolid - tVoid;
                                solidConnectionFlux[connection.id] = -t * deltaT;
                            }
                        }
                    }

                    solidConductionSource[scv.dofIndex()] +=
                        couplingManager->conductionSource(CouplingManager::solidDomainIdx, element, solidFVGeometry, solidElemVolVars, scv)*scv.volume();

                    solidConvectionSource[scv.dofIndex()] +=
                        couplingManager->convectionSource(CouplingManager::solidDomainIdx, element, solidFVGeometry, solidElemVolVars, scv)*scv.volume();
                }
            }
        }

        std::vector<Scalar> voidEnergyOutflowFlux(voidGridGeometry->numDofs(), 0.0);

        for (const auto& element : elements(voidGridView))
        {
            voidFVGeometry.bind(element);
            voidElemVolVars.bind(element, voidFVGeometry, sol[CouplingManager::voidDomainIdx]);
            couplingManager->bindCouplingContext(CouplingManager::voidDomainIdx, element);

            for (const auto& scv : scvs(voidFVGeometry))
            {
                if (couplingManager->isCoupledPore(CouplingManager::voidDomainIdx, scv.dofIndex()))
                {
                    // iterate over all connection throats
                    const auto& connections = couplingManager->couplingMapper().voidToSolidConnections(scv.dofIndex());
                    const auto& contextForPore = couplingManager->couplingContext()[scv.dofIndex()];

                    for (const auto& [localConnectionIdx, connection] : enumerate(connections))
                    {
                        if (!voidConnectionVisited[connection.id])
                        {
                            voidConnectionVisited[connection.id] = true;
                            const auto& context = contextForPore[localConnectionIdx];

                            if (voidProblem->sourceMode() == decltype(voidProblem->sourceMode())::convection)
                                voidConnectionFlux[connection.id] = couplingManager->convectiveHeatFluxForOneConnection(connection, context);
                            else
                            {
                                const Scalar t = couplingManager->getConnectionTransmissiblity(CouplingManager::voidDomainIdx, connection, voidElemVolVars, scv);
                                const Scalar deltaT = sol[CouplingManager::solidDomainIdx][connection.solidVertexIdx][0] - sol[CouplingManager::voidDomainIdx][connection.voidVertexIdx][1];
                                voidConnectionFlux[connection.id] = t * deltaT;
                            }
                        }
                    }

                    voidConductionSource[scv.dofIndex()] += couplingManager->conductionSource(CouplingManager::voidDomainIdx, element,
                                                                                              voidFVGeometry, voidElemVolVars, scv)
                                                                                              * scv.volume();

                    voidConvectionSource[scv.dofIndex()] += couplingManager->convectionSource(CouplingManager::voidDomainIdx, element,
                                                                                              voidFVGeometry, voidElemVolVars, scv)
                                                                                              * scv.volume();
                }

                if (voidProblem->onOutletBoundary(scv))
                    voidEnergyOutflowFlux[scv.dofIndex()] += voidProblem->heatOutFlowCondition(element, voidFVGeometry, voidElemVolVars, scv) * scv.volume();
            }
        }

        // write vtk output
        voidVtkWriter.write(1);
        solidVtkWriter.write(1);
        couplingWriter.write("coupling");

        static const bool fluxVerbose = getParam<bool>("Problem.BoundaryFluxVerbose", false);

        voidProblem->setOutletToDirichletForOutput(sol[CouplingManager::voidDomainIdx]);
        solidProblem->setInletToDirichletForOutput(sol[CouplingManager::solidDomainIdx]);

        const auto voidInletFlux = voidBoundaryFlux.getFlux(std::vector<int>{voidProblem->inletPoreLabel()}, fluxVerbose);
        const auto voidOutletFlux = voidBoundaryFlux.getFlux(std::vector<int>{voidProblem->outletPoreLabel()}, fluxVerbose);
        const auto voidHeaterFlux = voidBoundaryFlux.getFlux(std::vector<int>{voidProblem->heaterPoreLabel()}, fluxVerbose);
        const auto solidInletFlux = solidBoundaryFlux.getFlux(std::vector<int>{voidProblem->inletPoreLabel()}, fluxVerbose);
        const auto solidOutletFlux = solidBoundaryFlux.getFlux(std::vector<int>{voidProblem->outletPoreLabel()}, fluxVerbose);
        const auto solidHeaterFlux = solidBoundaryFlux.getFlux(std::vector<int>{voidProblem->heaterPoreLabel()}, fluxVerbose);
        using std::abs;

        std::cout << "\n\n ***** Boundary fluxes \n\n" << std::endl;
        std::cout << "Void inlet flux (total - mass/heat): " << voidInletFlux << std::endl;
        std::cout << "Void outlet flux (total - mass/heat): " << voidOutletFlux << std::endl;
        std::cout << "Delta heat flux abs(outlet) - abs(inlet): " << abs(voidOutletFlux.totalFlux[1]) - abs(voidInletFlux.totalFlux[1]) << "\n\n" << std::endl;

        std::cout << "Void heater flux (total - mass/heat): " << voidHeaterFlux << "\n\n" << std::endl;

        std::cout << "Solid inlet flux (total - heat): " << solidInletFlux << std::endl;
        std::cout << "Solid outlet flux (total - heat): " << solidOutletFlux << std::endl;
        std::cout << "Solid inlet Robin flux: " << totalSolidInletHeatFlux << std::endl;
        std::cout << "Delta heat flux abs(outlet) - abs(inlet): " << abs(solidOutletFlux.totalFlux[1]) - abs(solidInletFlux.totalFlux[1]) << "\n\n" << std::endl;

        std::cout << "Solid heater flux (total - heat): " << solidHeaterFlux << "\n\n" << std::endl;

        Scalar totalConvectiveHeatFluxIn = 0.0;
        Scalar totalConvectiveHeatFluxOut = 0.0;
        Scalar voidConductionInletHeatFlux = 0.0;
        std::vector<bool> visited(voidGridGeometry->numDofs(), false);

        for (const auto& element : elements(voidGridView))
        {
            voidFVGeometry.bind(element);
            voidElemVolVars.bind(element, voidFVGeometry, sol[CouplingManager::voidDomainIdx]);
            couplingManager->bindCouplingContext(CouplingManager::voidDomainIdx, element);

            for (const auto& scv : scvs(voidFVGeometry))
            {
                const auto dofIdx = scv.dofIndex();

                if (voidInletFlux.fluxPerPore.count(dofIdx))
                    voidConductionInletHeatFlux += voidProblem->robinInletConductiveHeatFlux(element, voidFVGeometry, voidElemVolVars, scv) * scv.volume();

                if (visited[dofIdx])
                    continue;
                else
                    visited[dofIdx] = true;

                if (voidInletFlux.fluxPerPore.count(dofIdx))
                {
                    const Scalar tIn = getParamFromGroup<Scalar>(voidProblem->paramGroup(), "Problem.InletTemperature");
                    const Scalar enthalypy = std::decay_t<decltype(voidElemVolVars)>::VolumeVariables::FluidSystem::enthalpy(tIn, voidElemVolVars[scv].pressure(0));
                    totalConvectiveHeatFluxIn += voidInletFlux.fluxPerPore.at(dofIdx)[0] * enthalypy;
                }
                if (voidOutletFlux.fluxPerPore.count(dofIdx))
                    totalConvectiveHeatFluxOut += voidOutletFlux.fluxPerPore.at(dofIdx)[0] * voidElemVolVars[scv].enthalpy(0);
            }
        }


        const Scalar sumCouplingVoid = std::accumulate(voidConnectionFlux.cbegin(), voidConnectionFlux.cend(), 0.0);
        const Scalar sumCouplingSolid = std::accumulate(solidConnectionFlux.cbegin(), solidConnectionFlux.cend(), 0.0);

        std::cout << "Sum of all void coupling heat fluxes: " << sumCouplingVoid << std::endl;
        std::cout << "Sum of all solid coupling heat fluxes: " << sumCouplingSolid << "\n\n" << std::endl;

        std::cout << "abs(Q_total,out,void) - abs(Q_total,in,void)- abs(Q_total,heater,void): " << abs(voidOutletFlux.totalFlux[1]) - abs(voidInletFlux.totalFlux[1]) -abs(voidHeaterFlux.totalFlux[1]) << std::endl;
        std::cout << "abs(Q_conv,out,void) - abs(Q_conv,in,void)- abs(Q_total,heater,void): " << abs(totalConvectiveHeatFluxOut) - abs(totalConvectiveHeatFluxIn) - abs(voidHeaterFlux.totalFlux[1]) << std::endl;


        std::cout << "Total sum void " << voidInletFlux[1] + voidOutletFlux[1] + voidHeaterFlux[1] - sumCouplingVoid << std::endl;

        std::cout << "\n\n ***Exchange *** " << std::endl;
        std::cout << "sumCouplingVoid: " << sumCouplingVoid << std::endl;
        std::cout << "solidInletFlux: " << solidInletFlux << std::endl;
        std::cout << "voidHeaterFlux: " << voidHeaterFlux[1] << std::endl;
        std::cout << "Total exchange flux (sumCouplingVoid + solidInletFlux - voidHeaterFlux) " << sumCouplingVoid + solidInletFlux[1] - voidHeaterFlux[1] << std::endl;

        std::cout << "\n\n *** Void convection *** " << std::endl;
        std::cout << "Void inlet convective heat flux (\u1E41 * h): " << totalConvectiveHeatFluxIn << std::endl;
        std::cout << "Void outlet convective heat flux (\u1E41 * h): " << totalConvectiveHeatFluxOut << std::endl;
        const Scalar sumConvection = totalConvectiveHeatFluxIn + totalConvectiveHeatFluxOut;
        std::cout << "abs((\u1E41 * h)_out) - abs((\u1E41 * h)_in): " << abs(totalConvectiveHeatFluxOut) - abs(totalConvectiveHeatFluxIn) << "\n\n" << std::endl;

        std::cout << "\n\n ***Inlet conduction ***" << std::endl;
        std::cout << "Void inlet conductive heat flux: " << -voidConductionInletHeatFlux << std::endl; // negative sign because value is used in source term
        std::cout << "Solid inlet Robin flux: " << -totalSolidInletHeatFlux << std::endl; // negative sign because value is used in source term
        const Scalar sumInletCond = -voidConductionInletHeatFlux - totalSolidInletHeatFlux;
        std::cout << "Sum " << std::abs(voidConductionInletHeatFlux) + std::abs(totalSolidInletHeatFlux) << std::endl;

        std::cout << "\n\n ***Heater conduction ***" << std::endl;
        std::cout << "Solid heater flux: " << solidHeaterFlux << std::endl;
        std::cout << "Void heater flux: " << voidHeaterFlux.totalFlux[1] << std::endl;
        const Scalar sumHeaterFlux = solidHeaterFlux[0] + voidHeaterFlux[1];
        std::cout << "Sum :" << solidHeaterFlux[0] + voidHeaterFlux[1] << std::endl;

        std::cout << "\n\n ***Sum of sums  *** " << std::endl;
        std::cout << sumConvection + sumInletCond + sumHeaterFlux << std::endl;

        Dune::VTKWriter fluxVtkWriter(voidGridGeometry->gridView());
        std::vector<Scalar> voidMassFlux(voidGridGeometry->numDofs());
        std::vector<Scalar> voidHeatFlux(voidGridGeometry->numDofs());
        std::vector<bool> isConsidered(voidGridGeometry->numDofs(), false);
        fluxVtkWriter.addVertexData(isConsidered, "isConsidered");
        fluxVtkWriter.addVertexData(voidMassFlux, "voidMassFlux");
        fluxVtkWriter.addVertexData(voidHeatFlux, "voidHeatFlux");
        fluxVtkWriter.addVertexData(voidEnergyOutflowFlux, "voidEnergyOutflowFlux");

        for (const auto& [dofIdx, value] : voidInletFlux.fluxPerPore)
        {
            isConsidered[dofIdx] = true;
            voidMassFlux[dofIdx] = value[0];
            voidHeatFlux[dofIdx] = value[1];
        }
        fluxVtkWriter.write("inlet_flux");

        std::fill(isConsidered.begin(), isConsidered.end(), false);
        std::fill(voidMassFlux.begin(), voidMassFlux.end(), 0.0);
        std::fill(voidHeatFlux.begin(), voidHeatFlux.end(), 0.0);
        for (const auto& [dofIdx, value] : voidOutletFlux.fluxPerPore)
        {
            isConsidered[dofIdx] = true;
            voidMassFlux[dofIdx] = value[0];
            voidHeatFlux[dofIdx] = value[1];
        }
        fluxVtkWriter.write("outlet_flux");
    }
    else
    {
        // time loop
        timeLoop->start(); do
        {
            // solve the non-linear system with time step control
            nonLinearSolver.solve(sol, *timeLoop);

            // make the new solution the old solution
            solOld = sol;
            voidGridVariables->advanceTimeStep();
            solidGridVariables->advanceTimeStep();

            // advance to the time loop to the next step
            timeLoop->advanceTimeStep();

            for (const auto& element : elements(voidGridView))
            {
                voidFVGeometry.bind(element);
                voidElemVolVars.bind(element, voidFVGeometry, sol[CouplingManager::voidDomainIdx]);
                couplingManager->bindCouplingContext(CouplingManager::voidDomainIdx, element);

                for (const auto& scv : scvs(voidFVGeometry))
                {
                    if (couplingManager->isCoupledPore(CouplingManager::voidDomainIdx, scv.dofIndex()))
                    {
                        voidConductionSource[scv.dofIndex()] += couplingManager->conductionSource(CouplingManager::voidDomainIdx, element,
                                                                                                  voidFVGeometry, voidElemVolVars, scv)
                                                                                                  * scv.volume();

                        voidConvectionSource[scv.dofIndex()] += couplingManager->convectionSource(CouplingManager::voidDomainIdx, element,
                                                                                                  voidFVGeometry, voidElemVolVars, scv)
                                                                                                  * scv.volume();
                    }
                }
            }

            for (const auto& element : elements(solidGridView))
            {
                solidFVGeometry.bind(element);
                solidElemVolVars.bind(element, solidFVGeometry, sol[CouplingManager::solidDomainIdx]);
                couplingManager->bindCouplingContext(CouplingManager::solidDomainIdx, element);

                for (const auto& scv : scvs(solidFVGeometry))
                {
                    if (couplingManager->isCoupledPore(CouplingManager::solidDomainIdx, scv.dofIndex()))
                    {
                        solidConductionSource[scv.dofIndex()] += couplingManager->conductionSource(CouplingManager::solidDomainIdx, element,
                                                                                                   solidFVGeometry, solidElemVolVars, scv)
                                                                                                   * scv.volume();

                        solidConvectionSource[scv.dofIndex()] += couplingManager->convectionSource(CouplingManager::solidDomainIdx, element,
                                                                                                   solidFVGeometry, solidElemVolVars, scv)
                                                                                                   * scv.volume();
                    }
                }
            }

            // write vtk output
            voidVtkWriter.write(timeLoop->time());
            solidVtkWriter.write(timeLoop->time());

            // report statistics of this time step
            timeLoop->reportTimeStep();

            // set new dt as suggested by newton solver
            timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

        } while (!timeLoop->finished());

        timeLoop->finalize(voidGridView.comm());
        timeLoop->finalize(solidGridView.comm());
    }

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
