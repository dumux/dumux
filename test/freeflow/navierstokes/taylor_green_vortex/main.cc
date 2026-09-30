// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Taylor-Green vortex test for the (hybrid) CVFE Navier-Stokes models
 *        (stationary and instationary, 2D and 3D), see README.md.
 */

#include <config.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/typetraits/griddiscretization.hh>

#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/extrusion.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>
#if SIMPLEX_GRID
#include <dumux/io/grid/gridmanager_alu.hh>
#endif
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>

// the new assembler relies on declarations of the old one
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/assembler.hh>

#include <dumux/freeflow/navierstokes/momentum/velocityoutput.hh>
#include <test/freeflow/navierstokes/analyticalsolutionvectors.hh>
#include <test/freeflow/navierstokes/errors_cvfe.hh>

#include "properties.hh"

namespace Dumux {

/*!
 * \brief Computes the numerical and the analytical kinetic energy
 *        \f$ E = \frac{1}{2} \int_\Omega \rho \|\mathbf{u}\|^2 \f$ by quadrature
 */
template<class Problem, class GridVariables, class SolutionVector>
std::pair<double, double> kineticEnergy(const Problem& problem,
                                        const GridVariables& gridVariables,
                                        const SolutionVector& x,
                                        double density,
                                        int order = 5)
{
    using GridGeometry = GridDiscretization_t<GridVariables>;
    using Extrusion = Extrusion_t<GridGeometry>;
    const auto& gridDiscretization = Dumux::gridDiscretization(problem);
    auto fvGeometry = localView(gridDiscretization);
    auto curGridVars = [&]() -> decltype(auto)
    {
        if constexpr (requires { gridVariables.curGridVolVars(); })
            return gridVariables.curGridVolVars();
        else
            return gridVariables.curGridVars();
    };
    auto elemVars = localView(curGridVars());

    double energy = 0.0;
    double exactEnergy = 0.0;
    for (const auto& element : elements(gridDiscretization.gridView()))
    {
        fvGeometry.bind(element);
        elemVars.bind(element, fvGeometry, x);
        const auto geometry = fvGeometry.elementGeometry();
        const auto elemSol = elementSolution(element, elemVars, fvGeometry);
        const auto& quad = Dune::QuadratureRules<double, GridGeometry::GridView::dimension>::rule(geometry.type(), order);
        for (auto&& qp : quad)
        {
            const auto weight = qp.weight() * Extrusion::integrationElement(geometry, qp.position());
            const auto velocity = evalSolutionAtLocalPos(element, geometry, gridDiscretization, elemSol, qp.position());
            const auto exactVelocity = problem.analyticalSolution(geometry.global(qp.position()));
            energy += 0.5*density*velocity.two_norm2()*weight;
            exactEnergy += 0.5*density*exactVelocity.two_norm2()*weight;
        }
    }

    return {energy, exactEnergy};
}

/*!
 * \brief Writes the errors and the kinetic energy of every time step to <Problem.Name>_errors.csv
 */
class TaylorGreenErrorWriter
{
public:
    explicit TaylorGreenErrorWriter(const std::string& name)
    : file_(name + "_errors.csv")
    {
        file_ << "t,dt,numDofsVelocity,numDofsPressure,h,"
              << "velocityL2,velocityH1,pressureL2,pressureH1,"
              << "kineticEnergy,kineticEnergyExact\n";
        file_ << std::scientific << std::setprecision(12);
    }

    template<class MomentumProblem, class MassProblem,
             class MomentumGridVariables, class MassGridVariables, class SolutionVector>
    void write(const MomentumProblem& momentumProblem,
               const MassProblem& massProblem,
               const MomentumGridVariables& momentumGridVariables,
               const MassGridVariables& massGridVariables,
               const SolutionVector& x,
               double time, double dt)
    {
        using namespace Dune::Indices;
        const auto& momentumX = x[_0];
        const auto& massX = x[_1];

        const auto [volume, velocityErrors] = calculateL2AndH1Errors(momentumProblem, momentumGridVariables, momentumX);
        const auto [massVolume, pressureErrors] = calculateL2AndH1Errors(massProblem, massGridVariables, massX);
        const auto density = getParam<double>("Component.LiquidDensity");
        const auto [energy, exactEnergy] = kineticEnergy(momentumProblem, momentumGridVariables, momentumX, density);

        const auto& gridDiscretization = Dumux::gridDiscretization(momentumProblem);
        static constexpr int dim = std::decay_t<decltype(gridDiscretization)>::GridView::dimension;
        const double h = std::pow(volume/gridDiscretization.gridView().size(0), 1.0/dim);

        std::cout << "[Errors] t = " << time
                  << " velocity L2 = " << velocityErrors[0] << " H1 = " << velocityErrors[1]
                  << " pressure L2 = " << pressureErrors[0] << " H1 = " << pressureErrors[1]
                  << " kinetic energy = " << energy << " (exact " << exactEnergy << ")" << std::endl;

        file_ << time << "," << dt << ","
              << gridDiscretization.numDofs() << "," << Dumux::gridDiscretization(massProblem).numDofs() << ","
              << h << ","
              << velocityErrors[0] << "," << velocityErrors[1] << ","
              << pressureErrors[0] << "," << pressureErrors[1] << ","
              << energy << "," << exactEnergy << std::endl;
    }

private:
    std::ofstream file_;
};

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tags for this problem
    using MomentumTypeTag = Properties::TTag::TYPETAG_MOMENTUM;
    using MassTypeTag = Properties::TTag::TaylorGreenTestMassBox;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // create the grid
    GridManager<GetPropType<MomentumTypeTag, Properties::Grid>> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometries
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    // the coupling manager
    using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the time loop (with a constant time step size)
    using Scalar = GetPropType<MomentumTypeTag, Properties::Scalar>;
    const bool isStationary = getParam<bool>("Problem.IsStationary");
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(getParam<Scalar>("TimeLoop.MaxTimeStepSize", dt));

    // the problems (initial and boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);
    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    // the solution vector
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    momentumProblem->applyInitialSolution(x[momentumIdx]);
    massProblem->applyInitialSolution(x[massIdx]);
    auto xOld = x;

    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    if (isStationary)
        couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    else
        couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x, xOld);

    massGridVariables->init(x[massIdx]);
    momentumGridVariables->init(x[momentumIdx]);

    // initialize the vtk output module
    const bool enableVtkOutput = getParam<bool>("Problem.EnableVtkOutput", true);
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    NavierStokesTest::AnalyticalSolutionVectors analyticalSolVectors(momentumProblem, massProblem);
    vtkWriter.addField(analyticalSolVectors.analyticalPressureSolution(), "pressureExact");
    vtkWriter.addField(analyticalSolVectors.analyticalVelocitySolution(), "velocityExact");
    if (enableVtkOutput)
        vtkWriter.write(0.0);

    // the errors (and kinetic energy) of every time step
    const bool printErrors = getParam<bool>("Problem.PrintErrors", true);
    TaylorGreenErrorWriter errorWriter(massProblem->name());
    const auto writeErrors = [&](Scalar time, Scalar timeStepSize)
    {
        if (printErrors)
            errorWriter.write(*momentumProblem, *massProblem, *momentumGridVariables, *massGridVariables,
                              x, time, timeStepSize);
    };

    // the assembler
    using Assembler = Experimental::MultiDomainAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = isStationary ?
        std::make_shared<Assembler>(
            std::make_tuple(momentumProblem, massProblem),
            std::make_tuple(momentumGridGeometry, massGridGeometry),
            std::make_tuple(momentumGridVariables, massGridVariables),
            couplingManager
        )
        :
        std::make_shared<Assembler>(
            std::make_tuple(momentumProblem, massProblem),
            std::make_tuple(momentumGridGeometry, massGridGeometry),
            std::make_tuple(momentumGridVariables, massGridVariables),
            couplingManager, timeLoop, xOld
        );

    // the linear solver
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    Dune::Timer timer;
    if (isStationary)
    {
        nonLinearSolver->solve(x);
        writeErrors(0.0, 0.0);

        analyticalSolVectors.update();
        if (enableVtkOutput)
            vtkWriter.write(1.0);
    }
    else
    {
        writeErrors(0.0, dt);

        timeLoop->start(); do
        {
            const Scalar newTime = timeLoop->time() + timeLoop->timeStepSize();

            // set the correct time level for the problem's boundary conditions
            momentumProblem->updateTime(newTime);
            massProblem->updateTime(newTime);

            // solve the non-linear system with time step control
            nonLinearSolver->solve(x, *timeLoop);
            xOld = x;

            // make the new solution the old solution
            momentumGridVariables->advanceTimeStep();
            massGridVariables->advanceTimeStep();

            writeErrors(newTime, timeLoop->timeStepSize());

            // advance the time loop to the next step
            timeLoop->advanceTimeStep();
            analyticalSolVectors.update(timeLoop->time());

            // write vtk output
            if (enableVtkOutput)
                vtkWriter.write(timeLoop->time());

            // report statistics of this time step
            timeLoop->reportTimeStep();

            // keep the time step size constant (for the temporal convergence study)
            timeLoop->setTimeStepSize(dt);

        } while (!timeLoop->finished());

        timeLoop->finalize(leafGridView.comm());
    }

    timer.stop();
    std::cout << "Simulation took " << timer.elapsed() << " seconds." << std::endl;

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
