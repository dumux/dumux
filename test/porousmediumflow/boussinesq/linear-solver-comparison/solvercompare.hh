// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Shared helper for benchmarking ISTL linear solvers on an assembled
 *        Boussinesq Jacobian (either formulation).
 *
 * For a fixed matrix/rhs pair (taken from a real Newton assembly during the
 * transient runs in main_generate_{pressure,vorticity}.cc and read back in by
 * main_solve_{pressure,vorticity}.cc), each candidate solver is asked to
 * solve A*x = b from x0 = 0 exactly once. This isolates the linear-algebra
 * question ("which solver/preconditioner handles this matrix well?") from
 * Newton/time-stepping behaviour, which would otherwise confound the
 * comparison.
 */
#ifndef DUMUX_BOUSSINESQ_SOLVERCOMPARE_HH
#define DUMUX_BOUSSINESQ_SOLVERCOMPARE_HH

#include <iostream>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

namespace Dumux {

struct SolverBenchmarkResult
{
    std::string state;
    std::string solver;
    bool converged = false;
    int iterations = -1;
    double reduction = 0.0;
    double solverElapsed = 0.0; //!< time reported by the solver itself (setup+solve)
    double wallTime = 0.0;      //!< time measured around the solve() call
    std::string note;
};

/*!
 * \brief Solve A*x = b once with LinearSolver, starting from x = 0.
 *
 * Never throws: solver setup/convergence failures (e.g. AMG on an
 * indefinite matrix, ILU hitting a zero pivot) are caught and reported
 * as a non-converged result instead of aborting the whole comparison.
 */
template<class LinearSolver, class Matrix, class Vector>
SolverBenchmarkResult benchmarkOneSolver(const std::string& state,
                                         const std::string& solverName,
                                         const Matrix& A,
                                         const Vector& b)
{
    SolverBenchmarkResult r;
    r.state = state;
    r.solver = solverName;

    try
    {
        Matrix ACopy(A);
        Vector bCopy(b);
        Vector x(b.size());
        x = 0.0;

        LinearSolver solver;
        Dune::Timer timer;
        const auto result = solver.solve(ACopy, x, bCopy);
        r.wallTime = timer.elapsed();

        r.converged = static_cast<bool>(result.converged);
        r.iterations = result.iterations;
        r.reduction = result.reduction;
        r.solverElapsed = result.elapsed;
    }
    catch (const Dune::Exception& e)
    {
        r.converged = false;
        r.note = e.what();
    }
    catch (const std::exception& e)
    {
        r.converged = false;
        r.note = e.what();
    }

    return r;
}

/*!
 * \brief Run the full candidate-solver set against one assembled system.
 *
 * \tparam Assembler the FVAssembler instantiation the Jacobian/residual came from
 */
template<class Assembler>
std::vector<SolverBenchmarkResult>
benchmarkAllSolvers(const std::string& state,
                    typename Assembler::JacobianMatrix& A,
                    typename Assembler::ResidualType& b)
{
    using LSTraits = SeqLinearSolverTraits;
    using LATraits = LinearAlgebraTraitsFromAssembler<Assembler>;

    std::vector<SolverBenchmarkResult> results;

    results.push_back(benchmarkOneSolver<UMFPackIstlSolver<LSTraits, LATraits>>(
        state, "UMFPack (direct)", A, b));
    results.push_back(benchmarkOneSolver<AMGBiCGSTABIstlSolver<LSTraits, LATraits>>(
        state, "AMG-BiCGSTAB", A, b));
    results.push_back(benchmarkOneSolver<AMGCGIstlSolver<LSTraits, LATraits>>(
        state, "AMG-CG", A, b));
    results.push_back(benchmarkOneSolver<ILUBiCGSTABIstlSolver<LSTraits, LATraits>>(
        state, "ILU-BiCGSTAB", A, b));
    results.push_back(benchmarkOneSolver<ILURestartedGMResIstlSolver<LSTraits, LATraits>>(
        state, "ILU-GMRes", A, b));
    results.push_back(benchmarkOneSolver<SSORBiCGSTABIstlSolver<LSTraits, LATraits>>(
        state, "SSOR-BiCGSTAB", A, b));

    return results;
}

inline void printSolverComparison(const std::string& title,
                                  const std::vector<SolverBenchmarkResult>& results)
{
    std::cout << "\n=== " << title << " ===\n";
    std::cout << std::left  << std::setw(18) << "solver"
              << std::left  << std::setw(12) << "state"
              << std::right << std::setw(10) << "converged"
              << std::right << std::setw(12) << "iterations"
              << std::right << std::setw(14) << "reduction"
              << std::right << std::setw(16) << "solver time [s]"
              << std::right << std::setw(14) << "wall time [s]" << "\n";

    for (const auto& r : results)
    {
        std::cout << std::left  << std::setw(18) << r.solver
                  << std::left  << std::setw(12) << r.state
                  << std::right << std::setw(10) << (r.converged ? "yes" : "NO")
                  << std::right << std::setw(12) << r.iterations
                  << std::scientific << std::setprecision(3)
                  << std::right << std::setw(14) << r.reduction
                  << std::fixed << std::setprecision(4)
                  << std::right << std::setw(16) << r.solverElapsed
                  << std::right << std::setw(14) << r.wallTime;
        if (!r.note.empty())
            std::cout << "  (" << r.note << ")";
        std::cout << "\n";
    }
}

inline void writeSolverComparisonCsvHeader(std::ostream& os)
{
    os << "formulation,state,solver,converged,iterations,reduction,solver_time_s,wall_time_s,note\n";
}

inline void writeSolverComparisonCsv(std::ostream& os,
                                     const std::string& formulation,
                                     const std::vector<SolverBenchmarkResult>& results)
{
    os << std::scientific << std::setprecision(6);
    for (const auto& r : results)
    {
        os << formulation << "," << r.state << "," << r.solver << ","
           << (r.converged ? 1 : 0) << "," << r.iterations << ","
           << r.reduction << "," << r.solverElapsed << "," << r.wallTime << ","
           << "\"" << r.note << "\"\n";
    }
}

} // end namespace Dumux

#endif
