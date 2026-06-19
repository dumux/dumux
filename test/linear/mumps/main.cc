// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Standalone test for the direct MUMPS solver backend (sequential path).
// Builds a known 1D Poisson (tridiagonal) system A x = b and checks the solution.
#include <config.h>

#include <iostream>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/mumpssolver.hh>

#if DUMUX_HAVE_MUMPS
using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
struct LATraits { using Matrix = ::Matrix; using Vector = ::Vector; };
#endif

int main(int argc, char* argv[])
{
    using namespace Dumux;
    Dune::MPIHelper::instance(argc, argv);

#if DUMUX_HAVE_MUMPS
    constexpr std::size_t n = 100;

    // Tridiagonal 1D Poisson matrix: diag 2, off-diag -1.
    Matrix A;
    A.setBuildMode(Matrix::random);
    A.setSize(n, n);
    for (std::size_t i = 0; i < n; ++i)
        A.setrowsize(i, (i == 0 || i == n - 1) ? 2 : 3);
    A.endrowsizes();
    for (std::size_t i = 0; i < n; ++i)
    {
        if (i > 0)     A.addindex(i, i - 1);
                       A.addindex(i, i);
        if (i < n - 1) A.addindex(i, i + 1);
    }
    A.endindices();
    for (std::size_t i = 0; i < n; ++i)
    {
        if (i > 0)     A[i][i - 1] = -1.0;
                       A[i][i]     =  2.0;
        if (i < n - 1) A[i][i + 1] = -1.0;
    }

    // Pick a known solution xExact, compute b = A xExact.
    Vector xExact(n), b(n), x(n);
    for (std::size_t i = 0; i < n; ++i)
        xExact[i] = std::sin(0.1 * static_cast<double>(i)) + 1.0;
    A.mv(xExact, b);
    x = 0.0;

    using Solver = DirectSolverMumps<SeqLinearSolverTraits, LATraits>;

    Solver solver;

    // Solve twice to exercise the analysis-reuse path (second solve reuses symbolic factorization).
    solver.solve(A, x, b);
    x = 0.0;
    solver.solve(A, x, b);

    double maxErr = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        maxErr = std::max(maxErr, std::abs(x[i] - xExact[i]));

    std::cout << "MUMPS sequential solve: max error = " << maxErr << std::endl;
    if (maxErr > 1e-10)
        DUNE_THROW(Dune::Exception, "MUMPS solution incorrect, max error " << maxErr);

    std::cout << "Test passed." << std::endl;
    return 0;
#else
    std::cout << "MUMPS not available, skipping test." << std::endl;
    return 77; // CMake/dune skip code
#endif
}
