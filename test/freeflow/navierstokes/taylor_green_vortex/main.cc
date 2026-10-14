// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Taylor-Green vortex test for the (hybrid) CVFE Navier-Stokes models
 *        (stationary and instationary, 2D and 3D), see README.md.
 */

#include <config.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/typetraits/griddiscretization.hh>

#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/method.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>
#if SIMPLEX_GRID
#include <dumux/io/grid/gridmanager_alu.hh>
#endif
#include <dumux/io/vtkoutputmodule.hh>
#if DUMUX_HAVE_GRIDFORMAT
#include <dumux/io/gridwriter.hh>
#include <dumux/io/cvfegridfunction.hh>
#endif

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>
// the new assembler relies on declarations of the old one
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/assembler.hh>
#include <dumux/multidomain/multistagemultidomainassembler.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>

#include <dumux/freeflow/navierstokes/momentum/velocityoutput.hh>
#include <test/freeflow/navierstokes/analyticalsolutionvectors.hh>
#include <test/freeflow/navierstokes/errors.hh>

#include "properties.hh"

namespace Dumux {

/*!
 * \brief The vorticity from the velocity gradients (of all primary variables):
 *        the scalar \f$ \omega = \partial_x v - \partial_y u \f$ in 2D and
 *        the vector \f$ \boldsymbol{\omega} = \nabla \times \mathbf{u} \f$ in 3D
 */
template<class Indices, class Gradients>
auto vorticityFromGradients(const Gradients& gradients)
{
    static constexpr int dim = Gradients::value_type::dimension;
    const auto du = [&](int i, int j) { return gradients[Indices::velocity(i)][j]; };
    if constexpr (dim == 2)
        return du(1, 0) - du(0, 1);
    else
        return Dune::FieldVector<double, 3>{ du(2, 1) - du(1, 2), du(0, 2) - du(2, 0), du(1, 0) - du(0, 1) };
}

//! The squared norm of a scalar or vector
template<class T>
double squaredNorm(const T& v)
{
    if constexpr (std::is_arithmetic_v<T>)
        return v*v;
    else
        return v.two_norm2();
}

//! The benchmark indicators based on the velocity
struct TaylorGreenIndicators
{
    double kineticEnergy = 0.0;
    double kineticEnergyExact = 0.0;
    double enstrophy = 0.0;
    double enstrophyExact = 0.0;
};

/*!
 * \brief Computes the numerical and the analytical kinetic energy
 *        \f$ E = \frac{1}{2} \int_\Omega \rho \|\mathbf{u}\|^2 \, \mathrm{d}x \f$, enstrophy
 *        \f$ \boldsymbol{\varepsilon} = \frac{1}{2} \int_\Omega \rho \|\boldsymbol{\omega}\|^2 \, \mathrm{d}x \f$
 *        by quadrature
 */
template<class Problem, class GridVariables, class SolutionVector>
TaylorGreenIndicators computeIndicators(const Problem& problem,
                                        const GridVariables& gridVariables,
                                        const SolutionVector& x,
                                        double density,
                                        int order)
{
    using GridGeometry = GridDiscretization_t<GridVariables>;
    using Extrusion = Extrusion_t<GridGeometry>;
    using Indices = typename Problem::Indices;
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

    TaylorGreenIndicators indicators;
    for (const auto& element : elements(gridDiscretization.gridView()))
    {
        fvGeometry.bind(element);
        elemVars.bind(element, fvGeometry, x);
        const auto geometry = fvGeometry.elementGeometry();
        const auto elemSol = elementSolution(element, elemVars, fvGeometry);
        const auto& quad = Dune::QuadratureRules<double, GridGeometry::GridView::dimension>::rule(geometry.type(), order);
        for (auto&& qp : quad)
        {
            const auto& localPos = qp.position();
            const auto globalPos = geometry.global(localPos);
            const auto weight = qp.weight() * Extrusion::integrationElement(geometry, localPos);

            const auto velocity = evalSolutionAtLocalPos(element, geometry, gridDiscretization, elemSol, localPos);
            const auto exactVelocity = problem.analyticalSolution(globalPos);
            indicators.kineticEnergy += 0.5*density*velocity.two_norm2()*weight;
            indicators.kineticEnergyExact += 0.5*density*exactVelocity.two_norm2()*weight;

            const auto gradients = evalGradientsAtLocalPos(element, geometry, gridDiscretization, elemSol, localPos);
            const auto vorticity = vorticityFromGradients<Indices>(gradients);
            const auto exactVorticity = problem.analyticalVorticity(globalPos);
            indicators.enstrophy += 0.5*density*squaredNorm(vorticity)*weight;
            indicators.enstrophyExact += 0.5*density*squaredNorm(exactVorticity)*weight;
        }
    }

    return indicators;
}

/*!
 * \brief The vorticity of the velocity of the momentum discretization as grid function
 *        (in the sense of dune-functions, for the higher-order output)
 * \note The vorticity is discontinuous across elements; at nodes shared by several elements,
 *       the output contains the value of one of the elements.
 */
template<class GridGeometry, class SolutionVector, class Indices>
class VorticityGridFunction
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    class LocalFunction
    {
    public:
        LocalFunction(const GridGeometry& gg, const SolutionVector& sol)
        : gridGeometry_(&gg), sol_(&sol)
        {}

        void bind(const Element& element)
        { element_.emplace(element); }

        void unbind()
        { element_.reset(); }

        auto operator()(const typename Element::Geometry::LocalCoordinate& localPos) const
        {
            const auto elemSol = elementSolution(*element_, *sol_, *gridGeometry_);
            const auto gradients = evalGradientsAtLocalPos(*element_, element_->geometry(), *gridGeometry_, elemSol, localPos);
            return vorticityFromGradients<Indices>(gradients);
        }

    private:
        const GridGeometry* gridGeometry_;
        const SolutionVector* sol_;
        std::optional<Element> element_;
    };

    VorticityGridFunction(const GridGeometry& gg, const SolutionVector& sol)
    : gridGeometry_(&gg), sol_(&sol)
    {}

    friend LocalFunction localFunction(const VorticityGridFunction& f)
    { return LocalFunction(*f.gridGeometry_, *f.sol_); }

private:
    const GridGeometry* gridGeometry_;
    const SolutionVector* sol_;
};

/*!
 * \brief Writes the errors, the kinetic energy and the enstrophy of every time step to a csv file
 *        (Problem.ErrorFile, default: <Problem.Name>_errors.csv)
 *
 * The mesh sizes of velocity and pressure
 * \f$ h^v = (|\Omega|/N_v)^{1/d} \f$ and \f$ h^p = (|\Omega|/N_p)^{1/d} \f$ are written, where
 * \f$ N_v \f$ and \f$ N_p \f$ are the numbers of velocity and pressure degrees of freedom
 * (as in Schneider & Koch \cite Schneider2024). They allow to compare schemes with different element types
 * and numbers of degrees of freedom per element.
 */
class TaylorGreenErrorWriter
{
public:
    explicit TaylorGreenErrorWriter(const std::string& fileName)
    : file_(fileName)
    {
        file_ << "t,dt,numDofsVelocity,numDofsPressure,hVelocity,hPressure,"
              << "velocityL2,velocityH1,pressureL2,pressureH1,"
              << "kineticEnergy,kineticEnergyExact,enstrophy,enstrophyExact\n";
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

        const auto [volume, velocityErrors] = NavierStokesTest::calculateL2AndH1Errors(momentumProblem, momentumGridVariables, momentumX, quadratureOrder_);
        const auto [massVolume, pressureErrors] = NavierStokesTest::calculateL2AndH1Errors(massProblem, massGridVariables, massX, quadratureOrder_);
        const auto density = getParam<double>("Component.LiquidDensity");
        const auto indicators = computeIndicators(momentumProblem, momentumGridVariables, momentumX, density, quadratureOrder_);

        const auto& gridDiscretization = Dumux::gridDiscretization(momentumProblem);
        const auto numDofsVelocity = gridDiscretization.numDofs();
        const auto numDofsPressure = Dumux::gridDiscretization(massProblem).numDofs();
        static constexpr int dim = std::decay_t<decltype(gridDiscretization)>::GridView::dimension;
        const double hVelocity = std::pow(volume/numDofsVelocity, 1.0/dim);
        const double hPressure = std::pow(massVolume/numDofsPressure, 1.0/dim);

        std::cout << "[Errors] t = " << time
                  << " velocity L2 = " << velocityErrors[0] << " H1 = " << velocityErrors[1]
                  << " pressure L2 = " << pressureErrors[0] << " H1 = " << pressureErrors[1] << "\n"
                  << "[Indicators] t = " << time
                  << " kinetic energy = " << indicators.kineticEnergy << " (exact " << indicators.kineticEnergyExact << ")"
                  << " enstrophy = " << indicators.enstrophy << " (exact " << indicators.enstrophyExact << ")" << std::endl;

        file_ << time << "," << dt << ","
              << numDofsVelocity << "," << numDofsPressure << ","
              << hVelocity << "," << hPressure << ","
              << velocityErrors[0] << "," << velocityErrors[1] << ","
              << pressureErrors[0] << "," << pressureErrors[1] << ","
              << indicators.kineticEnergy << "," << indicators.kineticEnergyExact << ","
              << indicators.enstrophy << "," << indicators.enstrophyExact << std::endl;
    }

private:
    // the quadrature error of the error norms is of the same order as the squared discretization
    // error itself for rules of order 5 (with a considerable constant for the higher-order schemes),
    // so a rule of order 2k + 2 for the polynomial degree k = 3 of hybrid PQ3 is used
    static constexpr int quadratureOrder_ = 8;
    std::ofstream file_;
};

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tags for this problem
    using MomentumTypeTag = Properties::TTag::TYPETAG_MOMENTUM;
    using MassTypeTag = Properties::TTag::TYPETAG_MASS;

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

    // the vorticity at the element centers (a scalar in 2D, a vector in 3D)
    using MomentumIndices = typename MomentumProblem::Indices;
    using Vorticity = typename MomentumProblem::Vorticity;
    std::vector<Vorticity> vorticity(leafGridView.size(0)), vorticityExact(leafGridView.size(0));
    vtkWriter.addField(vorticity, "vorticity", Vtk::FieldType::element);
    vtkWriter.addField(vorticityExact, "vorticityExact", Vtk::FieldType::element);
    const auto updateVorticity = [&](Scalar time)
    {
        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = momentumGridGeometry->elementMapper().index(element);
            const auto geometry = element.geometry();
            const auto& center = referenceElement(geometry).position(0, 0);
            const auto elemSol = elementSolution(element, x[momentumIdx], *momentumGridGeometry);
            const auto gradients = evalGradientsAtLocalPos(element, geometry, *momentumGridGeometry, elemSol, center);
            vorticity[eIdx] = vorticityFromGradients<MomentumIndices>(gradients);
            vorticityExact[eIdx] = momentumProblem->analyticalVorticity(geometry.global(center), time);
        }
    };

    // higher-order output on Lagrange cells of the polynomial degree of the velocity, which shows
    // the actual velocity field of the momentum discretization (the velocity of the vtk output
    // above is evaluated at the element centers)
#if DUMUX_HAVE_GRIDFORMAT
    static constexpr int outputOrder = MomentumGridGeometry::discMethod == DiscretizationMethods::pq3 ? 3 : 2;
    using HigherOrderWriter = IO::GridWriter<std::decay_t<decltype(leafGridView)>, outputOrder>;
    using VorticityFunction = VorticityGridFunction<MomentumGridGeometry, std::decay_t<decltype(x[momentumIdx])>, MomentumIndices>;
    std::unique_ptr<HigherOrderWriter> higherOrderWriter;
    if (enableVtkOutput)
    {
        const auto format = IO::Format::pvd_with(IO::Format::vtu.with({
            .encoder = IO::Encoding::base64,
            .compressor = IO::Compression::none,
            .data_format = IO::VTK::DataFormat::inlined
        }));
        higherOrderWriter = std::make_unique<HigherOrderWriter>(
            format, leafGridView, massProblem->name() + "_ho", IO::order<outputOrder>
        );
        higherOrderWriter->setPointField("velocity", IO::cvfeGridFunction(*momentumGridGeometry, x[momentumIdx]));
        higherOrderWriter->setPointField("pressure", IO::cvfeGridFunction(*massGridGeometry, x[massIdx]));
        higherOrderWriter->setPointField("vorticity", VorticityFunction(*momentumGridGeometry, x[momentumIdx]));
    }
#endif

    const auto writeOutput = [&](Scalar time)
    {
        if (!enableVtkOutput)
            return;

        updateVorticity(time);
        vtkWriter.write(time);
#if DUMUX_HAVE_GRIDFORMAT
        higherOrderWriter->write(time);
#endif
    };
    writeOutput(0.0);

    // the errors and the benchmark indicators of every time step
    const bool printErrors = getParam<bool>("Problem.PrintErrors", true);
    // the error file can be written to another directory than the ParaView output
    TaylorGreenErrorWriter errorWriter(getParam<std::string>("Problem.ErrorFile", massProblem->name() + "_errors.csv"));
    const auto writeErrors = [&](Scalar time, Scalar timeStepSize)
    {
        if (printErrors)
            errorWriter.write(*momentumProblem, *massProblem, *momentumGridVariables, *massGridVariables,
                              x, time, timeStepSize);
    };

    Dune::Timer timer;
    if (isStationary)
    {
        // the assembler, linear and non-linear solver
        using Assembler = Experimental::MultiDomainAssembler<Traits, CouplingManager, DiffMethod::numeric>;
        auto assembler = std::make_shared<Assembler>(
            std::make_tuple(momentumProblem, massProblem),
            std::make_tuple(momentumGridGeometry, massGridGeometry),
            std::make_tuple(momentumGridVariables, massGridVariables),
            couplingManager
        );
        using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
        auto linearSolver = std::make_shared<LinearSolver>();
        using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
        auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

        nonLinearSolver->solve(x);
        writeErrors(0.0, 0.0);

        analyticalSolVectors.update();
        writeOutput(1.0);
    }
    else
    {
        // the time integration method (TimeLoop.Scheme), see test/freeflow/navierstokes/sincos
        const auto timeSteppingScheme = getParam<std::string>("TimeLoop.Scheme", "ImplicitEuler");
        std::shared_ptr<Experimental::MultiStageMethod<Scalar>> timeSteppingMethod;
        if (timeSteppingScheme == "ImplicitEuler")
            timeSteppingMethod = std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>();
        else if (timeSteppingScheme == "CrankNicolson")
            timeSteppingMethod = std::make_shared<Experimental::MultiStage::Theta<Scalar>>(0.5);
        else if (timeSteppingScheme == "DIRK2")
            timeSteppingMethod = std::make_shared<Experimental::MultiStage::DIRKSecondOrderAlexander<Scalar>>();
        else if (timeSteppingScheme == "DIRK3")
            timeSteppingMethod = std::make_shared<Experimental::MultiStage::DIRKThirdOrderAlexander<Scalar>>();
        else
            DUNE_THROW(ParameterException, "Unknown TimeLoop.Scheme " << timeSteppingScheme
                        << ". Use ImplicitEuler, CrankNicolson, DIRK2 or DIRK3.");

        // the assembler, linear and non-linear solver; the multi-stage multi-domain assembler does not
        // implement the residual assembly (assembleResidual), so the Newton solver only uses the shift criterion
        using Assembler = Experimental::MultiStageMultiDomainAssembler<Traits, CouplingManager, DiffMethod::numeric>;
        auto assembler = std::make_shared<Assembler>(
            std::make_tuple(momentumProblem, massProblem),
            std::make_tuple(momentumGridGeometry, massGridGeometry),
            std::make_tuple(momentumGridVariables, massGridVariables),
            couplingManager, timeSteppingMethod, xOld
        );
        assembler->setLinearSystem();
        using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
        auto linearSolver = std::make_shared<LinearSolver>();
        using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
        auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);
        Experimental::MultiStageTimeStepper<NewtonSolver> timeStepper(nonLinearSolver, timeSteppingMethod);

        writeErrors(0.0, dt);

        timeLoop->start(); do
        {
            const Scalar newTime = timeLoop->time() + timeLoop->timeStepSize();

            // the assembler sets the problem time of each stage (and refers to xOld as previous solution)
            timeStepper.step(x, timeLoop->time(), timeLoop->timeStepSize());
            xOld = x;

            // set the new time level for the error evaluation and the output
            momentumProblem->updateTime(newTime);
            massProblem->updateTime(newTime);

            // make the new solution the old solution
            momentumGridVariables->advanceTimeStep();
            massGridVariables->advanceTimeStep();

            writeErrors(newTime, timeLoop->timeStepSize());

            // advance the time loop to the next step
            timeLoop->advanceTimeStep();
            analyticalSolVectors.update(timeLoop->time());

            // write vtk output
            writeOutput(timeLoop->time());

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
