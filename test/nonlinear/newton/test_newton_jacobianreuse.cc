//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test of Newton's method with reuse of the Jacobian (Newton.EnableJacobianReuse)
 */
#include <config.h>

#include <cmath>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "test_newton_common.hh"

namespace Dumux {

// records the assemblies of Jacobian and residual (J) and of the residual alone (R)
class CountingScalarAssembler : public MockScalarAssembler
{
public:
    using MockScalarAssembler::MockScalarAssembler;

    void assembleResidual(const Scalar& x)
    {
        if (failNextResidualAssembly)
        {
            failNextResidualAssembly = false;
            DUNE_THROW(NumericalProblem, "Residual assembly failed on request");
        }

        MockScalarAssembler::assembleResidual(x);
        ++numResidualAssemblies;
        assemblies += 'R';
    }

    void assembleJacobianAndResidual(const Scalar& x)
    {
        MockScalarAssembler::assembleJacobianAndResidual(x);
        ++numJacobianAssemblies;
        assemblies += 'J';
    }

    std::size_t numIterations() const
    { return numJacobianAssemblies + numResidualAssemblies; }

    std::size_t numJacobianAssemblies = 0;
    std::size_t numResidualAssemblies = 0;
    std::string assemblies;
    bool failNextResidualAssembly = false;
};

// an assembler of a time-integration stage, whose Jacobian depends on the weights of the storage and spatial terms
class StageScalarAssembler : public CountingScalarAssembler
{
public:
    using CountingScalarAssembler::CountingScalarAssembler;

    std::pair<Scalar, Scalar> currentStageWeights() const
    { return weights; }

    std::pair<Scalar, Scalar> weights = {1.0, 1.0};
};

// keeps the matrix given to setMatrix for the solves with solve(x, b)
class ReusingScalarLinearSolver : public MockScalarLinearSolver
{
public:
    using MockScalarLinearSolver::solve;

    void setMatrix(const double& A)
    {
        A_ = A;
        ++numSetMatrix;
    }

    bool solve(double& x, const double& b) const
    {
        x = b/A_;
        return true;
    }

    std::size_t numSetMatrix = 0;

private:
    double A_ = 0.0;
};

template<class Assembler, class LinearSolver = ReusingScalarLinearSolver>
using Solver = NewtonSolver<Assembler, LinearSolver, DefaultPartialReassembler>;

void check(bool condition, const std::string& message)
{
    if (!condition)
        DUNE_THROW(Dune::Exception, message);
}

// x^2 - 5 = 0 from far away, where the Jacobian at the initial guess contracts poorly
void testSameRootAsNewton()
{
    const auto solve = [](const std::string& paramGroup)
    {
        auto assembler = std::make_shared<CountingScalarAssembler>();
        auto linearSolver = std::make_shared<ReusingScalarLinearSolver>();
        Solver<CountingScalarAssembler> solver(assembler, linearSolver, Dune::MPIHelper::getCommunication(), paramGroup);
        double x = 10.0;
        solver.solve(x);
        return std::make_tuple(x, assembler, linearSolver);
    };

    const auto [xNewton, assemblerNewton, linearSolverNewton] = solve("");
    const auto [xReuse, assemblerReuse, linearSolverReuse] = solve("Reuse");

    using std::abs; using std::sqrt;
    const double exact = sqrt(5.0);
    std::cout << std::setprecision(15)
              << "Newton: x = " << xNewton << ", " << assemblerNewton->numJacobianAssemblies << " Jacobians in "
              << assemblerNewton->numIterations() << " iterations\n"
              << "Reuse:  x = " << xReuse << ", " << assemblerReuse->numJacobianAssemblies << " Jacobians in "
              << assemblerReuse->numIterations() << " iterations (" << assemblerReuse->assemblies << ")" << std::endl;

    check(abs(xNewton - exact) < 1e-13*exact, "Newton's method did not find the root");
    check(abs(xReuse - exact) < 1e-8*exact, "Newton's method with Jacobian reuse did not find the root");
    check(assemblerNewton->numResidualAssemblies == 0, "Newton's method assembled the residual alone");
    check(linearSolverNewton->numSetMatrix == 0, "Newton's method set a matrix in the linear solver");
    check(assemblerReuse->numResidualAssemblies > 0, "No iteration reused the Jacobian");
    check(assemblerReuse->assemblies.starts_with("JJ"),
          "The Jacobian of the first iteration was reused before a contraction was observed: " + assemblerReuse->assemblies);
    check(linearSolverReuse->numSetMatrix == assemblerReuse->numJacobianAssemblies,
          "Not every assembled Jacobian was set once in the linear solver");
}

// a linear equation 2x - 6 = 0, whose Jacobian does not depend on the solution
template<class Assembler>
std::shared_ptr<Assembler> makeLinearAssembler()
{ return std::make_shared<Assembler>([](double x){ return 2.0*x - 6.0; }, [](double){ return 2.0; }); }

void testReuseAcrossSolves()
{
    auto assembler = makeLinearAssembler<StageScalarAssembler>();
    Solver<StageScalarAssembler> solver(assembler, std::make_shared<ReusingScalarLinearSolver>(),
                                        Dune::MPIHelper::getCommunication(), "Reuse");
    const auto checkAssemblies = [&](const std::string& expected, const std::string& situation)
    {
        check(assembler->assemblies == expected,
              situation + ": assemblies " + assembler->assemblies + " instead of " + expected);
    };

    double x = 0.0;
    for (int i = 0; i < 3; ++i)
        solver.solve(x);
    check(x == 3.0, "Wrong solution of the linear equation");
    checkAssemblies("JJRRRR", "Equal stage weights");

    assembler->weights = {1.0, 0.5};
    solver.solve(x);
    checkAssemblies("JJRRRRJJ", "Changed stage weights");

    assembler->failNextResidualAssembly = true;
    bool failed = false;
    try { solver.solve(x); }
    catch (const NumericalProblem&) { failed = true; }
    check(failed, "The requested failure of the residual assembly did not fail the solve");

    solver.solve(x);
    checkAssemblies("JJRRRRJJJJ", "Solve after a failed solve");
}

void testNoReuseAcrossSolvesWithoutStageWeights()
{
    auto assembler = std::make_shared<CountingScalarAssembler>();
    Solver<CountingScalarAssembler> solver(assembler, std::make_shared<ReusingScalarLinearSolver>(),
                                           Dune::MPIHelper::getCommunication(), "Reuse");

    double x = 10.0;
    solver.solve(x);
    const auto firstSolve = assembler->assemblies;
    check(firstSolve.find('R') != std::string::npos, "No reuse within a solve without stage weights: " + firstSolve);

    solver.solve(x);
    check(assembler->assemblies.substr(firstSolve.size()).starts_with("J"),
          "Without stage weights, a solve started with a Jacobian of the previous solve: " + assembler->assemblies);
}

template<class LinearSolver>
void testRejectedConfiguration(const std::string& paramGroup, const std::string& configuration)
{
    bool thrown = false;
    try
    {
        Solver<CountingScalarAssembler, LinearSolver> solver(std::make_shared<CountingScalarAssembler>(),
                                                             std::make_shared<LinearSolver>(),
                                                             Dune::MPIHelper::getCommunication(), paramGroup);
    }
    catch (const Dune::InvalidStateException&) { thrown = true; }
    check(thrown, "Accepted " + configuration);
}

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);

    Parameters::init([](Dune::ParameterTree& params){
        params["Reuse.Newton.EnableJacobianReuse"] = "true";
        params["PartialReassembly.Newton.EnableJacobianReuse"] = "true";
        params["PartialReassembly.Newton.EnablePartialReassembly"] = "true";
    });

    testSameRootAsNewton();
    testReuseAcrossSolves();
    testNoReuseAcrossSolvesWithoutStageWeights();
    testRejectedConfiguration<ReusingScalarLinearSolver>("PartialReassembly", "Jacobian reuse with partial reassembly");
    testRejectedConfiguration<MockScalarLinearSolver>("Reuse", "Jacobian reuse with a linear solver without setMatrix");

    return 0;
}
