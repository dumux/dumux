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
 * \ingroup NavierStokesTests
 * \brief Channel flow test for the staggered grid (Navier-)Stokes model.
 */

#include <config.h>

#include <ctime>
#include <iostream>
#include <fstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
// #include <dumux/freeflow/navierstokes/staggered/fluxoversurface.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/freeflow/navierstokes/velocityoutput.hh>

#include "problem_new.hh"
#include "../../analyticalsolution.hh"

#include <dumux/freeflow/navierstokes/newtonsolver.hh>

template<class Assembler, class SolutionVector, class LocalAssembler, class CouplingManager>
auto assembleBoundaryFluxes(const Assembler& assembler, const SolutionVector& curSol,
        CouplingManager& couplingManager)
{
    using MassSolutionVector = std::remove_reference_t<std::remove_cv_t<decltype(curSol[Dune::index_constant<1>()])>>;
    typename MassSolutionVector::block_type flux(0.0);

    for (const auto& element : elements(assembler.gridView(Dune::index_constant<1>())))
    {
        LocalAssembler localAssembler(assembler, element, curSol, couplingManager);
        localAssembler.bindLocalViews();

        for (const auto& scvf : scvfs(localAssembler.fvGeometry()))
        {
            if (scvf.boundary())
            {
                flux += localAssembler.localResidual().evalFlux(localAssembler.problem(),
                        element,
                        localAssembler.fvGeometry(),
                        localAssembler.curElemVolVars(),
                        localAssembler.elemFluxVarsCache(),
                        scvf
                        )
                    * scvf.area();
            }
        }
    }
    return flux;
}

namespace Dumux::Properties{

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::ChannelTest>
{
private:
    using Traits = MultiDomainTraits<TTag::ChannelTestMomentum, TTag::ChannelTestMass>;
public:
    using type = StaggeredFreeFlowCouplingManager<Traits>;
};

}

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::ChannelTestMomentum;
    using MassTypeTag = Properties::TTag::ChannelTestMass;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<MomentumTypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    momentumGridGeometry->update();

    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);
    massGridGeometry->update();

    // the coupling manager
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using CouplingManager = StaggeredFreeFlowCouplingManager<Traits>;

    auto couplingManager = std::make_shared<CouplingManager>();

    // the problem (boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);

    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    // get some time loop parameters
    using Scalar = GetPropType<MomentumTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = getParam<Scalar>("Restart.Time", 0);

    // the solution vector
    // the solution vector
    constexpr auto momentumIdx = Dune::index_constant<0>();
    constexpr auto massIdx = Dune::index_constant<1>();
    using SolutionVector = typename Traits::SolutionVector;
    //auto xPtr = std::make_shared<SolutionVector>();
    //auto xOldPtr = std::make_shared<SolutionVector>();
    SolutionVector x;// = *xPtr;
    //*xOldPtr = x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());
    momentumProblem->applyInitialSolution(x[momentumIdx]);
    massProblem->applyInitialSolution(x[massIdx]);
    auto xOld = x;

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    //if (getParam<Scalar>("Problem.InletVelocity") > 1e-6)
    //    timeLoop->setCheckPoint({200.0, 210.0});
    if (hasParam("Problem.OutputInterval"))
        timeLoop->setPeriodicCheckPoint(getParam<Scalar>("Problem.OutputInterval"));

    massProblem->setTimeLoop(timeLoop);
    momentumProblem->setTimeLoop(timeLoop);

    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);

    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x, xOld);
    momentumGridVariables->init(x[momentumIdx]);
    massGridVariables->init(x[massIdx]);

    // initialize the vtk output module
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());

    const bool isStationary = getParam<bool>("Problem.IsStationary", false);
    if (!isStationary)
        vtkWriter.write(0);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = isStationary ? std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                                std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                                std::make_tuple(momentumGridVariables, massGridVariables),
                                                                couplingManager)
                                  : std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                                std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                                std::make_tuple(momentumGridVariables, massGridVariables),
                                                                couplingManager,
                                                                timeLoop, xOld);
    // the linear solver
    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = PhasefieldMultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // // set up two surfaces over which fluxes are calculated
    // FluxOverSurface<GridVariables,
    //                 SolutionVector,
    //                 GetPropType<TypeTag, Properties::ModelTraits>,
    //                 GetPropType<TypeTag, Properties::LocalResidual>> flux(*gridVariables, x);
    // using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    // using Element = typename GridView::template Codim<0>::Entity;

    // using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    // const Scalar xMin = gridGeometry->bBoxMin()[0];
    // const Scalar xMax = gridGeometry->bBoxMax()[0];
    // const Scalar yMin = gridGeometry->bBoxMin()[1];
    // const Scalar yMax = gridGeometry->bBoxMax()[1];

    // // The first surface shall be placed at the middle of the channel.
    // // If we have an odd number of cells in x-direction, there would not be any cell faces
    // // at the position of the surface (which is required for the flux calculation).
    // // In this case, we add half a cell-width to the x-position in order to make sure that
    // // the cell faces lie on the surface. This assumes a regular cartesian grid.
    // const Scalar planePosMiddleX = xMin + 0.5*(xMax - xMin);
    // int numCellsX = getParam<std::vector<int>>("Grid.Cells")[0];

    // const unsigned int refinement = getParam<unsigned int>("Grid.Refinement", 0);

    // numCellsX *= (1<<refinement);

    // const Scalar offsetX = (numCellsX % 2 == 0) ? 0.0 : 0.5*((xMax - xMin) / numCellsX);

    // const auto p0middle = GlobalPosition{planePosMiddleX + offsetX, yMin};
    // const auto p1middle = GlobalPosition{planePosMiddleX + offsetX, yMax};
    // flux.addSurface("middle", p0middle, p1middle);

    // // The second surface is placed at the outlet of the channel.
    // const auto p0outlet = GlobalPosition{xMax, yMin};
    // const auto p1outlet = GlobalPosition{xMax, yMax};
    // flux.addSurface("outlet", p0outlet, p1outlet);

    using MassLocalAssembler = SubDomainCCLocalAssembler<massIdx, MassTypeTag, Assembler,
          DiffMethod::numeric, /*implict*/true>;
    std::ofstream fout_scalar;
    fout_scalar.open(getParam<std::string>("Problem.Name", "scalars") + ".txt");
    //auto boundaryFlux = assembleBoundaryFluxes<Assembler, SolutionVector,
    //     MassLocalAssembler, CouplingManager>(*assembler, xOld, *couplingManager);
    massProblem->writeScalars(xOld[massIdx],/* boundaryFlux,*/ fout_scalar);

    if (isStationary)
    {
        nonLinearSolver.solve(x);
        if (momentumProblem->hasAnalyticalSolution())
        {
            const auto exactPressure = getScalarAnalyticalSolution(*massProblem)[GetPropType<MassTypeTag, Properties::ModelTraits>::Indices::pressureIdx];
            const auto exactVelocity = getVelocityAnalyticalSolution(*momentumProblem);
            vtkWriter.addField(exactPressure, "pressureExact");
            vtkWriter.addField(exactVelocity, "velocityExact");
            vtkWriter.write(1);
        }
    }
    else
    {
        // time loop
        assembler->assembleResidual(x);
        //auto res = assembler->residual();
        //auto res_mom = res[momentumIdx];
        //auto res_mas = res[massIdx];
        //std::ofstream fout_res;
        //fout.open("res_mom.txt");
        //Dune::printvector(fout, res_mom, "", "");
        //fout_res.close();
        //fout_res.open("res_mas.txt");
        //Dune::printvector(fout_res, res_mas, "", "");
        timeLoop->start(); do
        {
            // solve the non-linear system with time step control
            nonLinearSolver.solve(x, *timeLoop);
            //auto A = assembler->jacobian();
            //auto A_mom = A[momentumIdx][momentumIdx];
            //auto A_mas = A[massIdx][massIdx];
            //Dune::writeMatrixToMatlab(A_mom, "a_mom_matrix.mat");
            //Dune::writeMatrixToMatlab(A_mom, "a_mas_matrix.mat");
            //auto res = assembler->residual();
            //auto res_mas = res[massIdx];
            //Dune::printvector(fout_res, res_mas, "", "");
            //auto res_u = res_mas[Indices];


            // make the new solution the old solution
            xOld = x;
            momentumGridVariables->advanceTimeStep();
            massGridVariables->advanceTimeStep();

            // advance to the time loop to the next step
            timeLoop->advanceTimeStep();

            // write vtk output
            vtkWriter.write(timeLoop->time());

            // write volume, surface and relative mass
            auto boundaryFlux = assembleBoundaryFluxes<Assembler, SolutionVector,
                 MassLocalAssembler, CouplingManager>(*assembler, xOld, *couplingManager);
            massProblem->writeScalars(xOld[massIdx], boundaryFlux, fout_scalar);

            // // calculate and print mass fluxes over the planes
            // flux.calculateMassOrMoleFluxes();
            // if(GetPropType<TypeTag, Properties::ModelTraits>::enableEnergyBalance())
            // {
            //     std::cout << "mass / energy flux at middle is: " << flux.netFlux("middle") << std::endl;
            //     std::cout << "mass / energy flux at outlet is: " << flux.netFlux("outlet") << std::endl;
            // }
            // else
            // {
            //     std::cout << "mass flux at middle is: " << flux.netFlux("middle") << std::endl;
            //     std::cout << "mass flux at outlet is: " << flux.netFlux("outlet") << std::endl;
            // }

            // // calculate and print volume fluxes over the planes
            // flux.calculateVolumeFluxes();
            // std::cout << "volume flux at middle is: " << flux.netFlux("middle")[0] << std::endl;
            // std::cout << "volume flux at outlet is: " << flux.netFlux("outlet")[0] << std::endl;

            // report statistics of this time step
            timeLoop->reportTimeStep();

            // set new dt as suggested by newton solver
            timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

        } while (!timeLoop->finished());

        timeLoop->finalize(leafGridView.comm());
        //fout_res.close();
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
