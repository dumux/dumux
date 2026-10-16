// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Test for the sequential MUMPS solver backend: one solver instance solves a sequence of
 *        systems whose values, sparsity pattern and size change, and reuses a factorization for
 *        several right-hand sides, one by one and in one solve, single-type (also with 4x4 blocks)
 *        and MultiType; several solver instances keep their own factorizations at the same time
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/mumpssolver.hh>

namespace Dumux::Test {

using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
using MTMatrix = Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<Matrix, Matrix>,
                                            Dune::MultiTypeBlockVector<Matrix, Matrix>>;
using MTVector = Dune::MultiTypeBlockVector<Vector, Vector>;

using BlockMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 4, 4>>;
using BlockVector = Dune::BlockVector<Dune::FieldVector<double, 4>>;

struct LinearAlgebraTraits { using Matrix = Test::Matrix; using Vector = Test::Vector; };
struct MTLinearAlgebraTraits { using Matrix = MTMatrix; using Vector = MTVector; };
struct BlockLinearAlgebraTraits { using Matrix = BlockMatrix; using Vector = BlockVector; };

/*!
 * \brief A rows x cols matrix with the entries (i, i + d) for the given diagonal offsets d,
 *        the value \a diagonal on the diagonal and unsymmetric off-diagonal values
 */
Matrix bandMatrix(std::size_t rows, std::size_t cols, const std::vector<int>& offsets,
                  double diagonal, double scale)
{
    Dune::MatrixIndexSet pattern(rows, cols);
    for (std::size_t i = 0; i < rows; ++i)
        for (int d : offsets)
            if (const long j = static_cast<long>(i) + d; j >= 0 && j < static_cast<long>(cols))
                pattern.add(i, j);

    Matrix A;
    pattern.exportIdx(A);
    for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
        for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
        {
            const auto i = rowIt.index();
            const auto j = colIt.index();
            *colIt = (i == j) ? diagonal : -scale*(i < j ? 0.3 : 0.2);
        }
    return A;
}

Vector exactSolution(std::size_t size, double shift)
{
    Vector x(size);
    for (std::size_t i = 0; i < size; ++i)
        x[i] = 1.0 + std::sin(0.1*static_cast<double>(i) + shift);
    return x;
}

template<class Solver, class M, class V>
void solveAndCheck(Solver& solver, const M& A, const V& xExact, const std::string& name)
{
    V b = xExact;
    A.mv(xExact, b);
    V x = xExact;
    x = 0.0;
    solver.solve(A, x, b);

    x -= xExact;
    const double error = x.infinity_norm();
    std::cout << name << ": max error " << error << std::endl;
    if (!(error < 1e-12))
        DUNE_THROW(Dune::Exception, name << ": wrong solution, max error " << error);
}

/*!
 * \brief Factorizes A with setMatrix, changes A afterwards and solves for two right-hand sides,
 *        which has to give the solutions of the unchanged A
 */
template<class Solver, class M, class V>
void solveWithReusedFactorization(Solver& solver, M A, const V& xExact, const std::string& name)
{
    V b = xExact;
    A.mv(xExact, b);
    solver.setMatrix(A);
    A *= 2.0;

    for (const double scale : {1.0, -3.0})
    {
        V scaledB = b;
        scaledB *= scale;
        V x = xExact;
        x = 0.0;
        solver.solve(x, scaledB);

        x *= 1.0/scale;
        x -= xExact;
        const double error = x.infinity_norm();
        std::cout << name << ", right-hand side scaled by " << scale << ": max error " << error << std::endl;
        if (!(error < 1e-12))
            DUNE_THROW(Dune::Exception, name << ": wrong solution with the reused factorization, max error " << error);
    }

    V wrongSize;
    bool detected = false;
    try { solver.solve(wrongSize, wrongSize); }
    catch (const Dune::InvalidStateException&) { detected = true; }
    if (!detected)
        DUNE_THROW(Dune::Exception, name << ": solving with vectors of the wrong size did not throw");
}

//! Several solver instances with their own factorizations, solved in turns, one of them destroyed
void solveWithSeveralInstances()
{
    using Solver = DirectSolverMumps<SeqLinearSolverTraits, LinearAlgebraTraits>;
    std::vector<std::unique_ptr<Solver>> solvers;
    std::vector<Matrix> matrices;
    for (int i = 0; i < 4; ++i)
    {
        solvers.push_back(std::make_unique<Solver>());
        matrices.push_back(bandMatrix(90, 90, {-1, 0, 1}, 2.0 + i, 1.0 + 0.5*i));
        solvers.back()->setMatrix(matrices.back());
    }
    solvers[1].reset();

    for (int round = 0; round < 2; ++round)
    {
        for (std::size_t i = 0; i < solvers.size(); ++i)
        {
            if (!solvers[i])
                continue;

            const auto xExact = exactSolution(90, 0.1*i + round);
            Vector b = xExact;
            matrices[i].mv(xExact, b);
            Vector x = xExact;
            x = 0.0;
            solvers[i]->solve(x, b);

            x -= xExact;
            const double error = x.infinity_norm();
            std::cout << "instance " << i << ", round " << round << ": max error " << error << std::endl;
            if (!(error < 1e-12))
                DUNE_THROW(Dune::Exception, "instance " << i << ": wrong solution, max error " << error);
        }
    }
}

MTMatrix coupledMatrix(std::size_t n0, std::size_t n1, const std::vector<int>& couplingOffsets, double scale)
{
    MTMatrix A;
    using namespace Dune::Indices;
    A[_0][_0] = bandMatrix(n0, n0, {-1, 0, 1}, 2.0*scale, scale);
    A[_1][_1] = bandMatrix(n1, n1, {-2, 0, 2}, 2.0*scale, scale);
    A[_0][_1] = bandMatrix(n0, n1, couplingOffsets, -0.25*scale, 0.5*scale);
    A[_1][_0] = bandMatrix(n1, n0, couplingOffsets, -0.25*scale, 0.5*scale);
    return A;
}

MTVector coupledExactSolution(std::size_t n0, std::size_t n1, double shift = 0.0)
{
    MTVector x;
    using namespace Dune::Indices;
    x[_0] = exactSolution(n0, shift);
    x[_1] = exactSolution(n1, 1.0 + 2.0*shift);
    return x;
}

//! An unsymmetric, diagonally dominant matrix of 4x4 blocks on three block diagonals
BlockMatrix blockBandMatrix(std::size_t rows)
{
    Dune::MatrixIndexSet pattern(rows, rows);
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = (i > 0 ? i - 1 : 0); j < std::min(rows, i + 2); ++j)
            pattern.add(i, j);

    BlockMatrix A;
    pattern.exportIdx(A);
    for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
        for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
            for (int a = 0; a < 4; ++a)
                for (int c = 0; c < 4; ++c)
                    (*colIt)[a][c] = (rowIt.index() == colIt.index() && a == c) ? 8.0
                        : 0.1*std::sin(1.0 + rowIt.index() + 2.0*colIt.index() + 3.0*a + 5.0*c);
    return A;
}

BlockVector blockExactSolution(std::size_t size, double shift)
{
    BlockVector x(size);
    for (std::size_t i = 0; i < size; ++i)
        for (int k = 0; k < 4; ++k)
            x[i][k] = 1.0 + std::sin(0.1*static_cast<double>(i) + 0.7*k + shift);
    return x;
}

/*!
 * \brief Solves for the columns A xExact[j] in one solve; every column has to equal the exact
 *        solution and the single right-hand side solve of the same column
 */
template<class Solver, class M, class V>
void solveSeveralRightHandSides(Solver& solver, const M& A, const std::vector<V>& xExact, const std::string& name)
{
    std::vector<V> b(xExact), x(xExact);
    for (std::size_t j = 0; j < xExact.size(); ++j)
    {
        A.mv(xExact[j], b[j]);
        x[j] = 0.0;
    }
    solver.setMatrix(A);
    solver.solve(x, b);

    double maxToExact = 0.0, maxToSingle = 0.0;
    for (std::size_t j = 0; j < xExact.size(); ++j)
    {
        V single = xExact[j];
        single = 0.0;
        solver.solve(single, b[j]);
        single -= x[j];
        V error = x[j];
        error -= xExact[j];
        maxToExact = std::max(maxToExact, error.infinity_norm());
        maxToSingle = std::max(maxToSingle, single.infinity_norm());
    }
    std::cout << name << ", " << xExact.size() << " right-hand sides: max error " << maxToExact
              << ", max difference to single solves " << maxToSingle << std::endl;
    if (!(maxToExact < 1e-12) || !(maxToSingle < 1e-12))
        DUNE_THROW(Dune::Exception, name << ": wrong solution of several right-hand sides");

    std::vector<V> tooFew(xExact.begin(), xExact.end() - 1);
    bool detected = false;
    try { solver.solve(tooFew, b); }
    catch (const Dune::InvalidStateException&) { detected = true; }
    if (!detected)
        DUNE_THROW(Dune::Exception, name << ": solving with fewer solutions than right-hand sides did not throw");

    std::vector<V> none;
    solver.solve(none, none);
}

} // end namespace Dumux::Test

int main(int argc, char* argv[])
{
    using namespace Dumux;
    using namespace Dumux::Test;
    Dune::MPIHelper::instance(argc, argv);

    DirectSolverMumps<SeqLinearSolverTraits, LinearAlgebraTraits> solver;
    {
        Vector x(10), b(10);
        b = 1.0;
        bool detected = false;
        try { solver.solve(x, b); }
        catch (const Dune::InvalidStateException&) { detected = true; }
        if (!detected)
            DUNE_THROW(Dune::Exception, "Solving with solve(x, b) before setMatrix did not throw");
    }
    solveAndCheck(solver, bandMatrix(100, 100, {-1, 0, 1}, 2.0, 1.0), exactSolution(100, 0.0), "first solve");
    solveAndCheck(solver, bandMatrix(100, 100, {-1, 0, 1}, 2.0, 1.0), exactSolution(100, 0.5), "same matrix");
    solveAndCheck(solver, bandMatrix(100, 100, {-1, 0, 1}, 6.0, 3.0), exactSolution(100, 0.0), "new values");
    solveAndCheck(solver, bandMatrix(100, 100, {-2, 0, 2}, 2.0, 1.0), exactSolution(100, 0.0), "other pattern, same number of entries");
    solveAndCheck(solver, bandMatrix(100, 100, {-3, -1, 0, 1, 3}, 3.0, 1.0), exactSolution(100, 0.0), "more entries");
    solveAndCheck(solver, bandMatrix(150, 150, {-1, 0, 1}, 2.0, 1.0), exactSolution(150, 0.0), "larger system");
    solveAndCheck(solver, bandMatrix(80, 80, {-1, 0, 1}, 2.0, 1.0), exactSolution(80, 0.0), "smaller system");
    solveWithReusedFactorization(solver, bandMatrix(120, 120, {-2, -1, 0, 1, 2}, 4.0, 1.0), exactSolution(120, 0.3), "reused factorization");
    solveAndCheck(solver, bandMatrix(120, 120, {-1, 0, 1}, 2.0, 1.0), exactSolution(120, 0.0), "solve after the reused factorization");

    DirectSolverMumps<SeqLinearSolverTraits, MTLinearAlgebraTraits> mtSolver;
    solveAndCheck(mtSolver, coupledMatrix(60, 40, {0}, 1.0), coupledExactSolution(60, 40), "multitype: first solve");
    solveAndCheck(mtSolver, coupledMatrix(60, 40, {0}, 2.0), coupledExactSolution(60, 40), "multitype: new values");
    solveAndCheck(mtSolver, coupledMatrix(60, 40, {1}, 1.0), coupledExactSolution(60, 40), "multitype: other coupling pattern");
    solveAndCheck(mtSolver, coupledMatrix(30, 70, {0}, 1.0), coupledExactSolution(30, 70), "multitype: other subdomain sizes, same total size");
    solveAndCheck(mtSolver, coupledMatrix(50, 90, {0}, 1.0), coupledExactSolution(50, 90), "multitype: larger system");
    solveWithReusedFactorization(mtSolver, coupledMatrix(40, 60, {1}, 1.5), coupledExactSolution(40, 60), "multitype: reused factorization");

    {
        std::vector<Vector> columns;
        for (int j = 0; j < 5; ++j)
            columns.push_back(exactSolution(100, 0.4*j));
        solveSeveralRightHandSides(solver, bandMatrix(100, 100, {-2, -1, 0, 1, 3}, 4.0, 1.0), columns, "several right-hand sides");

        std::vector<MTVector> mtColumns;
        for (int j = 0; j < 4; ++j)
            mtColumns.push_back(coupledExactSolution(40, 60, 0.3*j));
        solveSeveralRightHandSides(mtSolver, coupledMatrix(40, 60, {0}, 1.0), mtColumns, "multitype: several right-hand sides");

        DirectSolverMumps<SeqLinearSolverTraits, BlockLinearAlgebraTraits> blockSolver;
        std::vector<BlockVector> blockColumns;
        for (int j = 0; j < 12; ++j)
            blockColumns.push_back(blockExactSolution(70, 0.25*j));
        solveSeveralRightHandSides(blockSolver, blockBandMatrix(70), blockColumns, "4x4 blocks: several right-hand sides");
    }

    solveWithSeveralInstances();

    return 0;
}
