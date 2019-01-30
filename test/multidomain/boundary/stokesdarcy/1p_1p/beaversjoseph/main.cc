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
 * \brief A test problem for the coupled Stokes/Darcy problem (1p)
 */
#include <config.h>

#include <ctime>
#include <iostream>
#include <fstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/istl/io.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/partial.hh>
#include <dumux/common/geometry/diameter.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/method.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

#include <dumux/multidomain/staggeredtraits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/multidomain/boundary/stokesdarcy/couplingmanager.hh>

#include <dumux/freeflow/navierstokes/staggered/fluxoversurface.hh>

#include "darcyproblem.hh"
#include "stokesproblem.hh"

namespace Dumux {
namespace Properties {

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::StokesOneP>
{
    using Traits = StaggeredMultiDomainTraits<TypeTag, TypeTag, Properties::TTag::DarcyOneP>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DarcyOneP>
{
    using Traits = StaggeredMultiDomainTraits<Properties::TTag::StokesOneP, Properties::TTag::StokesOneP, TypeTag>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

} // end namespace Properties
} // end namespace Dumux

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Define the sub problem type tags
    using StokesTypeTag = Properties::TTag::StokesOneP;
    using DarcyTypeTag = Properties::TTag::DarcyOneP;

    // use dune-subgrid to create the individual grids
    static constexpr int dim = 2;
    using HostGrid = Dune::YaspGrid<2, Dune::TensorProductCoordinates<double, dim> >;
    using HostGridManager = Dumux::GridManager<HostGrid>;
    HostGridManager hostGridManager;
    hostGridManager.init();
    auto& hostGrid = hostGridManager.grid();

    using GlobalPosition = Dune::FieldVector<double, dim>;

    const auto lowerLeftPM = getParamFromGroup<GlobalPosition>("Darcy", "Grid.LowerLeftPM");
    const auto upperRightPM = getParamFromGroup<GlobalPosition>("Darcy", "Grid.UpperRightPM");

    auto elementSelectorStokes = [&](const auto& element)
    {
        const auto x = element.geometry().center()[0];
        const auto y = element.geometry().center()[1];
        return x < lowerLeftPM[0] || x > upperRightPM[0] || y > upperRightPM[1];
    };

    auto elementSelectorDarcy = [&](const auto& element)
    {
        const auto x = element.geometry().center()[0];
        const auto y = element.geometry().center()[1];
        return x > lowerLeftPM[0] && x < upperRightPM[0] && y < upperRightPM[1];
    };

    // subgrid Pointer
    auto stokesGridPtr = SubgridGridCreator<HostGrid>::makeGrid(hostGrid, elementSelectorStokes, "Stokes");
    auto darcyGridPtr = SubgridGridCreator<HostGrid>::makeGrid(hostGrid, elementSelectorDarcy, "Darcy");

    // we compute on the leaf grid view
    const auto& darcyGridView = darcyGridPtr->leafGridView();
    const auto& stokesGridView = stokesGridPtr->leafGridView();

    // create the finite volume grid geometry
    using StokesFVGridGeometry = GetPropType<StokesTypeTag, Properties::FVGridGeometry>;
    auto stokesFvGridGeometry = std::make_shared<StokesFVGridGeometry>(stokesGridView);
    stokesFvGridGeometry->update();
    using DarcyFVGridGeometry = GetPropType<DarcyTypeTag, Properties::FVGridGeometry>;
    auto darcyFvGridGeometry = std::make_shared<DarcyFVGridGeometry>(darcyGridView);
    darcyFvGridGeometry->update();

    using Traits = StaggeredMultiDomainTraits<StokesTypeTag, StokesTypeTag, DarcyTypeTag>;

    // the coupling manager
    using CouplingManager = StokesDarcyCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>(stokesFvGridGeometry, darcyFvGridGeometry);

    // the indices
    constexpr auto stokesCellCenterIdx = CouplingManager::stokesCellCenterIdx;
    constexpr auto stokesFaceIdx = CouplingManager::stokesFaceIdx;
    constexpr auto darcyIdx = CouplingManager::darcyIdx;

    // get some time loop parameters
    using Scalar = GetPropType<StokesTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // the problem (initial and boundary conditions)
    using StokesProblem = GetPropType<StokesTypeTag, Properties::Problem>;
    auto stokesProblem = std::make_shared<StokesProblem>(stokesFvGridGeometry, couplingManager);
    using DarcyProblem = GetPropType<DarcyTypeTag, Properties::Problem>;
    auto darcyProblem = std::make_shared<DarcyProblem>(darcyFvGridGeometry, couplingManager);

    // the solution vector
    Traits::SolutionVector sol;
    sol[stokesCellCenterIdx].resize(stokesFvGridGeometry->numCellCenterDofs());
    sol[stokesFaceIdx].resize(stokesFvGridGeometry->numFaceDofs());
    sol[darcyIdx].resize(darcyFvGridGeometry->numDofs());

    // get a solution vector storing references to the two Stokes solution vectors
    auto stokesSol = partial(sol, stokesCellCenterIdx, stokesFaceIdx);

    // apply initial solution for instationary problems
    stokesProblem->applyInitialSolution(stokesSol);
    darcyProblem->applyInitialSolution(sol[darcyIdx]);

    auto solOld = sol;

    couplingManager->init(stokesProblem, darcyProblem, sol);

    // the grid variables
    using StokesGridVariables = GetPropType<StokesTypeTag, Properties::GridVariables>;
    auto stokesGridVariables = std::make_shared<StokesGridVariables>(stokesProblem, stokesFvGridGeometry);
    stokesGridVariables->init(stokesSol);
    using DarcyGridVariables = GetPropType<DarcyTypeTag, Properties::GridVariables>;
    auto darcyGridVariables = std::make_shared<DarcyGridVariables>(darcyProblem, darcyFvGridGeometry);
    darcyGridVariables->init(sol[darcyIdx]);

    // intialize the vtk output module
    const auto stokesName = getParam<std::string>("Problem.Name") + "_" + stokesProblem->name();
    const auto darcyName = getParam<std::string>("Problem.Name") + "_" + darcyProblem->name();

    StaggeredVtkOutputModule<StokesGridVariables, std::decay_t<decltype(stokesSol)>> stokesVtkWriter(*stokesGridVariables, stokesSol, stokesName);
    GetPropType<StokesTypeTag, Properties::VtkOutputFields>::initOutputModule(stokesVtkWriter);

    stokesVtkWriter.addField(stokesProblem->getAnalyticalVelocityX(), "analyticalV_x");
    stokesVtkWriter.addVolumeVariable([](const auto& v){ return v.density() - 1.23012;}, "dRho");
    stokesVtkWriter.addVolumeVariable([](const auto& v){ return v.pressure() - 1e5;}, "dP");

    stokesVtkWriter.write(0.0);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;

#if INSTATIONARY
    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
    timeLoop->setPeriodicCheckPoint(10.0);

    // the assembler for a stationary problem
    auto assembler = std::make_shared<Assembler>(std::make_tuple(stokesProblem, stokesProblem, darcyProblem),
                                                 std::make_tuple(stokesFvGridGeometry->cellCenterFVGridGeometryPtr(),
                                                                 stokesFvGridGeometry->faceFVGridGeometryPtr(),
                                                                 darcyFvGridGeometry),
                                                 std::make_tuple(stokesGridVariables->cellCenterGridVariablesPtr(),
                                                                 stokesGridVariables->faceGridVariablesPtr(),
                                                                 darcyGridVariables),
                                                 couplingManager, timeLoop);
#else
    auto assembler = std::make_shared<Assembler>(std::make_tuple(stokesProblem, stokesProblem, darcyProblem),
                                                 std::make_tuple(stokesFvGridGeometry->cellCenterFVGridGeometryPtr(),
                                                                 stokesFvGridGeometry->faceFVGridGeometryPtr(),
                                                                 darcyFvGridGeometry),
                                                 std::make_tuple(stokesGridVariables->cellCenterGridVariablesPtr(),
                                                                 stokesGridVariables->faceGridVariablesPtr(),
                                                                 darcyGridVariables),
                                                 couplingManager);
#endif

    VtkOutputModule<DarcyGridVariables, GetPropType<DarcyTypeTag, Properties::SolutionVector>> darcyVtkWriter(*darcyGridVariables, sol[darcyIdx], darcyName);
    using DarcyVelocityOutput = GetPropType<DarcyTypeTag, Properties::VelocityOutput>;
    darcyVtkWriter.addVelocityOutput(std::make_shared<DarcyVelocityOutput>(*darcyGridVariables));
    GetPropType<DarcyTypeTag, Properties::VtkOutputFields>::initOutputModule(darcyVtkWriter);
    darcyVtkWriter.addVolumeVariable([](const auto& v){ return v.density() - 1.23012;}, "dRho");
    darcyVtkWriter.addVolumeVariable([](const auto& v){ return v.pressure() - 1e5;}, "dP");
    darcyVtkWriter.write(0.0);

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // set up two surfaces over which fluxes are calculated
    FluxOverSurface<StokesGridVariables,
                    decltype(stokesSol),
                    GetPropType<StokesTypeTag, Properties::ModelTraits>,
                    GetPropType<StokesTypeTag, Properties::LocalResidual>> flux(*stokesGridVariables, stokesSol);

    const auto p0frontChannel = GlobalPosition{darcyFvGridGeometry->bBoxMin()[0], darcyFvGridGeometry->bBoxMax()[1]};;
    const auto p1frontChannel = GlobalPosition{darcyFvGridGeometry->bBoxMin()[0], stokesFvGridGeometry->bBoxMax()[1]};
    flux.addSurface("frontChannel", p0frontChannel, p1frontChannel);

    const Scalar xMiddleChannel = getParam<Scalar>("FluxOverSurface.PosX");
    const auto p0MiddleChannel = GlobalPosition{xMiddleChannel,  darcyFvGridGeometry->bBoxMax()[1]};
    const auto p1MiddleChannel = GlobalPosition{xMiddleChannel, stokesFvGridGeometry->bBoxMax()[1]};
    flux.addSurface("middleChannel", p0MiddleChannel, p1MiddleChannel);

    const auto p0backChannel = GlobalPosition{darcyFvGridGeometry->bBoxMax()[0],  darcyFvGridGeometry->bBoxMax()[1]};
    const auto p1backChannel = GlobalPosition{darcyFvGridGeometry->bBoxMax()[0], stokesFvGridGeometry->bBoxMax()[1]};
    flux.addSurface("backChannel", p0backChannel, p1backChannel);

    const auto p0Inlet = stokesFvGridGeometry->bBoxMin();
    const auto p1Inlet = GlobalPosition{stokesFvGridGeometry->bBoxMin()[0], stokesFvGridGeometry->bBoxMax()[1]};
    flux.addSurface("inlet", p0Inlet, p1Inlet);

    const auto p0outlet = GlobalPosition{stokesFvGridGeometry->bBoxMax()[0], stokesFvGridGeometry->bBoxMin()[1]};
    const auto p1outlet = stokesFvGridGeometry->bBoxMax();
    flux.addSurface("outlet", p0outlet, p1outlet);

    auto printAllFluxes = [&]()
    {
        using Indices = typename GetPropType<StokesTypeTag, Properties::ModelTraits>::Indices;

        Scalar fluxFront = 0.0;
        Scalar fluxTop = 0.0;
        Scalar fluxBack = 0.0;

        auto fvGeometry = localView(*stokesFvGridGeometry);
        auto elemVolVars = localView(stokesGridVariables->curGridVolVars());
        auto elemFaceVars = localView(stokesGridVariables->curGridFaceVars());
        for (auto&& element : elements(stokesFvGridGeometry->gridView()))
        {
            couplingManager->bindCouplingContext(stokesCellCenterIdx, element, *assembler);

            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, stokesSol);
            elemFaceVars.bind(element, fvGeometry, stokesSol);

            for (const auto scvf : scvfs(fvGeometry))
            {
                if (!couplingManager->isCoupledEntity(CouplingManager::stokesIdx, scvf))
                    continue;

                if (scvf.directionIndex() == 0 && scvf.directionSign() > 0)
                {
                    fluxFront += stokesProblem->neumann(element, fvGeometry,
                                                        elemVolVars, elemFaceVars,
                                                        scvf)[Indices::conti0EqIdx] * scvf.area();
                }

                if (scvf.directionIndex() == 1 && scvf.directionSign() < 0)
                {
                    fluxTop += stokesProblem->neumann(element, fvGeometry,
                                                      elemVolVars, elemFaceVars,
                                                      scvf)[Indices::conti0EqIdx] * scvf.area();
                }

                if (scvf.directionIndex() == 0 && scvf.directionSign() < 0)
                {
                    fluxBack += stokesProblem->neumann(element, fvGeometry,
                                                       elemVolVars, elemFaceVars,
                                                       scvf)[Indices::conti0EqIdx] * scvf.area();
                }
            }
        }

        flux.calculateMassOrMoleFluxes();
        std::cout << "***********************" << std::endl;
        std::cout << "mass fluxes channel:" << std::endl;
        std::cout << "inlet         = " << flux.netFlux("inlet")[0] << std::endl;
        std::cout << "outlet        = " << flux.netFlux("outlet")[0] << std::endl;
        std::cout << "frontChannel  = " << flux.netFlux("frontChannel")[0] << std::endl;
        std::cout << "middleChannel = " << flux.netFlux("middleChannel")[0] << std::endl;
        std::cout << "backChannel   = " << flux.netFlux("backChannel")[0] << std::endl;

        std::cout << "inlet + frontChannel + front pm   = " << flux.netFlux("inlet")[0]+flux.netFlux("frontChannel")[0]+fluxFront << std::endl;
        std::cout << "ratio q_channel/q_in =" <<  flux.netFlux("middleChannel")[0] / flux.netFlux("inlet")[0] << std::endl;
        std::cout << "***********************" << std::endl;

        std::cout << "***********************" << std::endl;
        std::cout << "mass fluxes pm " << std::endl;
        std::cout << "front = " << fluxFront << std::endl;
        std::cout << "top   = " << fluxTop << std::endl;
        std::cout << "back  = " << fluxBack << std::endl;
        std::cout << "sum   = " << fluxFront+fluxTop+fluxBack << std::endl;
        std::cout << "***********************" << std::endl;

        std::cout << "Re middle channel = " << flux.netFlux("middleChannel")[0] / (1.2 * 1.5e-5) << std::endl;

        std::ofstream file;
        std::string outname = getParam<std::string>("Problem.Outputname");
        outname.append("_summary.out");
        file.open(outname, std::ios::out | std::ios::app);
        if (file.fail())
            throw std::ios_base::failure(std::strerror(errno));

        double theta = getParam<Scalar>("Darcy.SpatialParams.PermeabilityAngle");
#if INSTATIONARY
        file << timeLoop->time() << "\t";
#endif
        file << theta  << "\t " << fluxFront  << "\t " << fluxTop  << "\t " << fluxBack << "\t " << flux.netFlux("middleChannel")[0] << std::endl;;
        file.close();


    };

#if INSTATIONARY

        // time loop
        timeLoop->start(); do
        {
            // set previous solution for storage evaluations
            assembler->setPreviousSolution(solOld);

            // solve the non-linear system with time step control
            nonLinearSolver.solve(sol, *timeLoop);

            // make the new solution the old solution
            solOld = sol;
            stokesGridVariables->advanceTimeStep();
            darcyGridVariables->advanceTimeStep();

            // advance to the time loop to the next step
            timeLoop->advanceTimeStep();

            // write vtk output
            if (timeLoop->isCheckPoint() || timeLoop->finished())
            {
                stokesVtkWriter.write(timeLoop->time());
                darcyVtkWriter.write(timeLoop->time());
            }

            printAllFluxes();

            // report statistics of this time step
            timeLoop->reportTimeStep();

            // set new dt as suggested by newton solver
            timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

        } while (!timeLoop->finished());

        timeLoop->finalize(stokesGridView.comm());
        timeLoop->finalize(darcyGridView.comm());

#else

    // solve the non-linear system
    nonLinearSolver.solve(sol);

    printAllFluxes();

    // write vtk output
    stokesVtkWriter.write(1.0);
    darcyVtkWriter.write(1.0);
#endif

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
} // end main
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
