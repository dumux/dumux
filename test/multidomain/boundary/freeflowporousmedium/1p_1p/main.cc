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
 * \ingroup BoundaryTests
 * \brief A test problem for the coupled Stokes/Darcy problem (1p).
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/method.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/freeflow/navierstokes/velocityoutput.hh>

#include "properties.hh"

int main(int argc, char** argv)
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
    using FreeFlowMomentumTypeTag = Properties::TTag::FreeFlowOnePMomentum;
    using FreeFlowMassTypeTag = Properties::TTag::FreeFlowOnePMass;
    using DarcyTypeTag = Properties::TTag::DarcyOneP;

    // try to create a grid (from the given grid file or the input file)
    // for both sub-domains
    using DarcyGridManager = Dumux::GridManager<GetPropType<DarcyTypeTag, Properties::Grid>>;
    DarcyGridManager darcyGridManager;
    darcyGridManager.init("Darcy"); // pass parameter group

    using FreeFlowGridManager = Dumux::GridManager<GetPropType<FreeFlowMomentumTypeTag, Properties::Grid>>;
    FreeFlowGridManager freeFlowGridManager;
    freeFlowGridManager.init("FreeFlow"); // pass parameter group

    // we compute on the leaf grid view
    const auto& darcyGridView = darcyGridManager.grid().leafGridView();
    const auto& freeFlowGridView = freeFlowGridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using FreeFlowMomentumFVGridGeometry = GetPropType<FreeFlowMomentumTypeTag, Properties::GridGeometry>;
    auto freeFlowMomentumFVGridGeometry = std::make_shared<FreeFlowMomentumFVGridGeometry>(freeFlowGridView);
    freeFlowMomentumFVGridGeometry->update();
    using FreeFlowMassFVGridGeometry = GetPropType<FreeFlowMassTypeTag, Properties::GridGeometry>;
    auto freeFlowMassFVGridGeometry = std::make_shared<FreeFlowMassFVGridGeometry>(freeFlowGridView);
    freeFlowMassFVGridGeometry->update();
    using DarcyFVGridGeometry = GetPropType<DarcyTypeTag, Properties::GridGeometry>;
    auto darcyFvGridGeometry = std::make_shared<DarcyFVGridGeometry>(darcyGridView);
    darcyFvGridGeometry->update();

    using Traits = MultiDomainTraits<FreeFlowMomentumTypeTag, FreeFlowMassTypeTag, DarcyTypeTag>;
    using CouplingManager = FreeFlowPorousMediumCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the indices
    constexpr auto freeFlowMomentumIndex = CouplingManager::freeFlowMomentumIndex;
    constexpr auto freeFlowMassIndex = CouplingManager::freeFlowMassIndex;
    constexpr auto porousMediumIndex = CouplingManager::porousMediumIndex;

    // the problem (initial and boundary conditions)
    using FreeFlowMomentumProblem = GetPropType<FreeFlowMomentumTypeTag, Properties::Problem>;
    auto freeFlowMomentumProblem = std::make_shared<FreeFlowMomentumProblem>(freeFlowMomentumFVGridGeometry, couplingManager);
    using FreeFlowMassProblem = GetPropType<FreeFlowMassTypeTag, Properties::Problem>;
    auto freeFlowMassProblem = std::make_shared<FreeFlowMassProblem>(freeFlowMassFVGridGeometry, couplingManager);
    using DarcyProblem = GetPropType<DarcyTypeTag, Properties::Problem>;
    auto darcyProblem = std::make_shared<DarcyProblem>(darcyFvGridGeometry, couplingManager);

    // the solution vector
    Traits::SolutionVector sol;
    sol[freeFlowMomentumIndex].resize(freeFlowMomentumFVGridGeometry->numDofs());
    sol[freeFlowMassIndex].resize(freeFlowMassFVGridGeometry->numDofs());
    sol[porousMediumIndex].resize(darcyFvGridGeometry->numDofs());

    // the grid variables
    using FreeFlowMomentumGridVariables = GetPropType<FreeFlowMomentumTypeTag, Properties::GridVariables>;
    auto freeFlowMomentumGridVariables = std::make_shared<FreeFlowMomentumGridVariables>(freeFlowMomentumProblem, freeFlowMomentumFVGridGeometry);
    using DarcyGridVariables = GetPropType<DarcyTypeTag, Properties::GridVariables>;
    using FreeFlowMassGridVariables = GetPropType<FreeFlowMassTypeTag, Properties::GridVariables>;
    auto freeFlowMassGridVariables = std::make_shared<FreeFlowMassGridVariables>(freeFlowMassProblem, freeFlowMassFVGridGeometry);
    using DarcyGridVariables = GetPropType<DarcyTypeTag, Properties::GridVariables>;
    auto darcyGridVariables = std::make_shared<DarcyGridVariables>(darcyProblem, darcyFvGridGeometry);

    couplingManager->init(freeFlowMomentumProblem, freeFlowMassProblem, darcyProblem,
                          std::make_tuple(freeFlowMomentumGridVariables, freeFlowMassGridVariables, darcyGridVariables),
                          sol);

    freeFlowMomentumGridVariables->init(sol[freeFlowMomentumIndex]);
    freeFlowMassGridVariables->init(sol[freeFlowMassIndex]);
    darcyGridVariables->init(sol[porousMediumIndex]);

    // initialize the vtk output module
    VtkOutputModule vtkWriterFF(*freeFlowMassGridVariables, sol[freeFlowMassIndex], freeFlowMassProblem->name());
    GetPropType<FreeFlowMassTypeTag, Properties::IOFields>::initOutputModule(vtkWriterFF); // Add model specific output fields
    vtkWriterFF.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<FreeFlowMassGridVariables>>());

    const bool addAnalyticSolutionForV = getParamFromGroup<bool>(freeFlowMassProblem->paramGroup(), "Problem.AddAnalyticSolutionForV", false);

    if (!freeFlowMassProblem->verticalFlow() && addAnalyticSolutionForV)
    {
        using Scalar = typename Traits::Scalar;
        using std::sqrt;
        std::vector<Scalar> analyticalVelocityX(freeFlowMassFVGridGeometry->gridView().size(0));
        const Scalar dPdX = -freeFlowMomentumProblem->pressureDifference() / (freeFlowMassFVGridGeometry->bBoxMax()[0] - freeFlowMassFVGridGeometry->bBoxMin()[0]);
        static const Scalar mu = getParam<Scalar>("Component.LiquidKinematicViscosity") * getParam<Scalar>("Component.LiquidDensity");
        static const Scalar alpha = getParam<Scalar>("Darcy.SpatialParams.AlphaBeaversJoseph");
        static const Scalar K = getParam<Scalar>("Darcy.SpatialParams.Permeability");
        static const Scalar sqrtK = sqrt(K);
        const Scalar sigma = (freeFlowMassFVGridGeometry->bBoxMax()[1] - freeFlowMassFVGridGeometry->bBoxMin()[1])/sqrtK;

        const Scalar uB =  -K/(2.0*mu) * ((sigma*sigma + 2.0*alpha*sigma) / (1.0 + alpha*sigma)) * dPdX;

        for (const auto& element : elements(freeFlowMassFVGridGeometry->gridView()))
        {
            const auto eIdx = freeFlowMassFVGridGeometry->gridView().indexSet().index(element);
            const Scalar y = element.geometry().center()[1] - freeFlowMassFVGridGeometry->bBoxMin()[1];

            const Scalar u = uB*(1.0 + alpha/sqrtK*y) + 1.0/(2.0*mu) * (y*y + 2*alpha*y*sqrtK) * dPdX;
            analyticalVelocityX[eIdx] = u;
        }

        vtkWriterFF.addField(analyticalVelocityX, "analyticalV_x");
    }

    VtkOutputModule pmVtkWriter(*darcyGridVariables, sol[porousMediumIndex],  darcyProblem->name());
    GetPropType<DarcyTypeTag, Properties::IOFields>::initOutputModule(pmVtkWriter);
    pmVtkWriter.addVelocityOutput(std::make_shared<GetPropType<DarcyTypeTag, Properties::VelocityOutput>>(*darcyGridVariables));

    // the assembler for a stationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(freeFlowMomentumProblem, freeFlowMassProblem, darcyProblem),
                                                 std::make_tuple(freeFlowMomentumFVGridGeometry,
                                                                 freeFlowMassFVGridGeometry,
                                                                 darcyFvGridGeometry),
                                                 std::make_tuple(freeFlowMomentumGridVariables,
                                                                 freeFlowMassGridVariables,
                                                                 darcyGridVariables),
                                                 couplingManager);

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    std::cout << "NumDofs Face " << freeFlowMomentumFVGridGeometry->numDofs() << std::endl;
    std::cout << "NumDofs CC " << freeFlowMassFVGridGeometry->numDofs() << std::endl;
    std::cout << "NumDofs Darcy " << darcyFvGridGeometry->numDofs() << std::endl;

    // solve the non-linear system
    nonLinearSolver.solve(sol);

    // write vtk output
    vtkWriterFF.write(1.0);
    pmVtkWriter.write(1.0);

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
