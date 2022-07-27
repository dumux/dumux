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
#include <config.h>

#include <iostream>
#include <random>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_ug.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/incompressiblestokessolver.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/freeflow/navierstokes/velocityoutput.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::DFGChannelTestMomentum;
    using MassTypeTag = Properties::TTag::DFGChannelTestMass;

    // initialize MPI, finalize is done automatically on exit
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // create a grid
    using Grid = GetPropType<MassTypeTag, Properties::Grid>;
    Dumux::GridManager<Grid> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    // the coupling manager
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using CouplingManager = StaggeredFreeFlowCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);
    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    // the solution vector
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());
    std::cout << "Total number of dofs: "
        << massGridGeometry->numDofs() + momentumGridGeometry->numDofs()*MomentumGridGeometry::GridView::dimension
        << std::endl;

    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    // compute coupling stencil and afterwards initialize grid variables (need coupling information)
    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    massGridVariables->init(x[massIdx]);
    momentumGridVariables->init(x[momentumIdx]);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager);

    // initialize the vtk output module
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    vtkWriter.write(0.0);

    // the linear solver
    using M = typename Assembler::JacobianMatrix; using V = typename Assembler::ResidualType;
    using LinearSolver = IncompressibleStokesSolver<M, V, MassGridGeometry, LinearSolverTraits<MomentumGridGeometry>>;
    //using LinearSolver = UMFPackBackend;
    V dirichletDofs;
    dirichletDofs[momentumIdx].resize(momentumGridGeometry->numDofs());
    dirichletDofs[massIdx].resize(massGridGeometry->numDofs());
    for (const auto& element : elements(leafGridView))
    {
        const auto fvGeometry = localView(*momentumGridGeometry).bind(element);
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (scvf.boundary())
            {
                const auto bcTypes = momentumProblem->boundaryTypes(element, scvf);
                for (int i = 0; i < bcTypes.size(); ++i)
                    if (bcTypes.isDirichlet(i))
                        dirichletDofs[momentumIdx][fvGeometry.scv(scvf.insideScvIdx()).dofIndex()][i] = 1.0;
            }
        }
    }

    auto linearSolver = std::make_shared<LinearSolver>(massGridGeometry, dirichletDofs);

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // linearize & solve
    Dune::Timer timer;
    nonLinearSolver.solve(x);

    // write vtk output
    vtkWriter.write(1.0);

    // update coupling manager for output
    couplingManager->updateSolution(x);

    // evaluate benchmark indicators (only valid with inertia term at Re=20)
    if (momentumProblem->enableInertiaTerms())
    {
        static constexpr double cDrafReference = 5.57953523384;
        static constexpr double cLiftReference = 0.010618948146;
        static constexpr double pDiffReference = 0.11752016697;
        const auto [cDrag, cLift] = momentumProblem->evalDragAndLiftCoefficient(*momentumGridVariables, x[momentumIdx]);
        std::cout << "cDrag: " << cDrag
                << " (reference: " << cDrafReference << ")"
                << "\n"
                << "cLift: " << cLift
                << " (reference: " << cLiftReference << ")"
                << std::endl;

        const auto pDiff = massProblem->evalPressureDifference(x[massIdx]);
        std::cout << "pDiff: " << pDiff
                << " (reference: " << pDiffReference << ")"
                << std::endl;

        if (getParam<bool>("Problem.CheckIndicators", true))
        {
            using std::abs;
            if (abs(cDrag - cDrafReference) > 0.002)
                DUNE_THROW(Dune::Exception,
                    "Deviation from reference too large for drag coefficient: "
                    << cDrag << " (ref: " << cDrafReference << ")"
                );
            if (abs(cLift - cLiftReference) > 0.002)
                DUNE_THROW(Dune::Exception,
                    "Deviation from reference too large for lift coefficient: "
                    << cLift << " (ref: " << cLiftReference << ")"
                );
            if (abs(pDiff - pDiffReference) > 0.003)
                DUNE_THROW(Dune::Exception,
                    "Deviation from reference too large for pressure difference: "
                    << pDiff << " (ref: " << pDiffReference << ")"
                );
        }
    }

    timer.stop();
    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

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
