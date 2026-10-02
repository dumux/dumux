// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Test for the sequential MUMPS solver backend: one solver instance solves a sequence of
 *        systems whose values, sparsity pattern and size change, single-type and MultiType
 */
#include <config.h>

#include <cmath>
#include <iostream>
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

struct LinearAlgebraTraits { using Matrix = Test::Matrix; using Vector = Test::Vector; };
struct MTLinearAlgebraTraits { using Matrix = MTMatrix; using Vector = MTVector; };

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

MTVector coupledExactSolution(std::size_t n0, std::size_t n1)
{
    MTVector x;
    using namespace Dune::Indices;
    x[_0] = exactSolution(n0, 0.0);
    x[_1] = exactSolution(n1, 1.0);
    return x;
}

} // end namespace Dumux::Test

int main(int argc, char* argv[])
{
    using namespace Dumux;
    using namespace Dumux::Test;
    Dune::MPIHelper::instance(argc, argv);

    DirectSolverMumps<SeqLinearSolverTraits, LinearAlgebraTraits> solver;
    solveAndCheck(solver, bandMatrix(100, 100, {-1, 0, 1}, 2.0, 1.0), exactSolution(100, 0.0), "first solve");
    solveAndCheck(solver, bandMatrix(100, 100, {-1, 0, 1}, 2.0, 1.0), exactSolution(100, 0.5), "same matrix");
    solveAndCheck(solver, bandMatrix(100, 100, {-1, 0, 1}, 6.0, 3.0), exactSolution(100, 0.0), "new values");
    solveAndCheck(solver, bandMatrix(100, 100, {-2, 0, 2}, 2.0, 1.0), exactSolution(100, 0.0), "other pattern, same number of entries");
    solveAndCheck(solver, bandMatrix(100, 100, {-3, -1, 0, 1, 3}, 3.0, 1.0), exactSolution(100, 0.0), "more entries");
    solveAndCheck(solver, bandMatrix(150, 150, {-1, 0, 1}, 2.0, 1.0), exactSolution(150, 0.0), "larger system");
    solveAndCheck(solver, bandMatrix(80, 80, {-1, 0, 1}, 2.0, 1.0), exactSolution(80, 0.0), "smaller system");

    DirectSolverMumps<SeqLinearSolverTraits, MTLinearAlgebraTraits> mtSolver;
    solveAndCheck(mtSolver, coupledMatrix(60, 40, {0}, 1.0), coupledExactSolution(60, 40), "multitype: first solve");
    solveAndCheck(mtSolver, coupledMatrix(60, 40, {0}, 2.0), coupledExactSolution(60, 40), "multitype: new values");
    solveAndCheck(mtSolver, coupledMatrix(60, 40, {1}, 1.0), coupledExactSolution(60, 40), "multitype: other coupling pattern");
    solveAndCheck(mtSolver, coupledMatrix(30, 70, {0}, 1.0), coupledExactSolution(30, 70), "multitype: other subdomain sizes, same total size");
    solveAndCheck(mtSolver, coupledMatrix(50, 90, {0}, 1.0), coupledExactSolution(50, 90), "multitype: larger system");

    return 0;
}
