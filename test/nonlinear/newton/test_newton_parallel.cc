//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <iostream>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <chrono>
#include <thread>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/istl/bvector.hh>

#include <dumux/common/timeloop.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "test_newton_common.hh"

/*
  Parallel Newton test. If the solver fails, we reduce the time step size.
  In parallel this can lead to dead locks if we don't correctly handle exceptions.
  We can't avoid all such situation if code in Dune throw and communicates after.
  However, we avoid it in situation where we can recover.
 */

namespace Dumux {

class MockScalarParallelLinearSolver
{
public:
    MockScalarLinearSolver(int rank, std::shared_ptr<TimeLoop<double>> timeLoop)
    : rank_(rank), timeLoop_(timeLoop) {}

    void setResidualReduction(double residualReduction) {}

    template<class Vector>
    bool solve(const double& A, Vector& x, const Vector& b) const
    {
        // solver construction
        auto solver = constructSolver_(A, b);

        // error handling: make sure the solver was successfully constructed on all processes
        // and throw on all processes if solver construction failed
        bool success = static_cast<bool>(solver);
        int successRemote = success;
        if (Dune::MPIHelper::instance().size() > 1)
            successRemote = Dune::MPIHelper::instance().getCommunication().min(success);

        if (!success)
            DUNE_THROW(Dune::Exception, "Could not create solver");
        else if (!successRemote)
            DUNE_THROW(Dune::Exception, "Could not create solver on remote process");

        // solver solve (here we assume that either all processes are successful or all fail)
        x = solver->solve();

        return true;
    }

    double norm(const double& residual) const
    {
        using std::abs;
        return abs(residual);
    }
private:
    int rank_;
    std::shared_ptr<TimeLoop<double>> timeLoop_;

    template<class Vector>
    struct Solver
    {
        Solver(const double& A, const Vector& b, int rank, const std::shared_ptr<TimeLoop<double>>& timeLoop)
        : A_(A), b_(b), rank_(rank), timeLoop_(timeLoop)
        {
            // constructor might fail and failure might be recoverable
            // this is what we are testing here
            if (rank_ == 1 && timeLoop_->timeStepSize() > 0.4 && timeLoop_->timeStepSize() < 0.9)
            {
                using namespace std::chrono_literals;
                std::this_thread::sleep_for(0.3s);
                DUNE_THROW(Dune::Exception, "This is a recoverable test error during solver constructor");
            }
        }

        auto solve() const
        {
            auto x = b_/A_;

            // collective communication emulating parallel solver
            x = Dune::MPIHelper::instance().getCommunication().min(x);

            // solver might not converge, this is recoverable and what we are testing here
            if (rank_ == 1 && timeLoop_->timeStepSize() > 0.9)
            {
                using namespace std::chrono_literals;
                std::this_thread::sleep_for(0.3s);
                DUNE_THROW(Dune::Exception, "This is a recoverable test error during solver solve");
            }

            return x;
        }
        const double& A_; const Vector& b_;
        int rank_; std::shared_ptr<TimeLoop<double>> timeLoop_;
    };

    template<class Vector>
    auto constructSolver_(const double& A, const Vector& b) const
    {
        try {
            return std::make_shared<Solver<Vector>>(A, b, rank_, timeLoop_);
        } catch (const Dune::Exception& e) {
            std::cerr << "Caught exception on solver construction: " << e.what() << std::endl;
            return std::decay_t<decltype(std::make_shared<Solver<Vector>>(A, b, rank_, timeLoop_))>();
        }
    }
};

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // initialize parameters
    Dumux::Parameters::init(argc, argv);

    // use the Newton solver to find a solution to a scalar equation
    using Assembler = MockScalarAssembler;
    using LinearSolver = MockScalarParallelLinearSolver;
    using Solver = NewtonSolver<Assembler, LinearSolver, DefaultPartialReassembler>;

    const int rank = Dune::MPIHelper::instance().rank();
    auto timeLoop = std::make_shared<TimeLoop<double>>(0.0, 1.0, 1.0);

    auto assembler = std::make_shared<Assembler>();
    auto linearSolver = std::make_shared<LinearSolver>(rank, timeLoop);
    auto solver = std::make_shared<Solver>(assembler, linearSolver);

    double initialGuess = 0.1;
    double x = initialGuess;

    if (rank == 0)
        std::cout << "Solving: x^2 - 5 = 0" << std::endl;

    solver->solve(x, *timeLoop);

    if (rank == 0)
    {
        // check if there have been two time step reductions
        if (Dune::FloatCmp::ne(timeLoop->timeStepSize(), 0.25, 1e-8))
            DUNE_THROW(Dune::Exception, "Time step reduction and recovery did not work");

        std::cout << "Solution: " << std::setprecision(15) << x
                << ", exact: " << std::sqrt(5.0)
                << ", error: " << std::abs(x-std::sqrt(5.0))/std::sqrt(5.0)*100 << "%" << std::endl;

        if (Dune::FloatCmp::ne(x, std::sqrt(5.0), 1e-13))
            DUNE_THROW(Dune::Exception, "Didn't find correct root: " << std::setprecision(15) << x << ", exact: " << std::sqrt(5.0));
    }

    return 0;
}
