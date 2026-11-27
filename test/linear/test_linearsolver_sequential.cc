//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief With SeqLinearSolverTraits every process solves its own system, also in a parallel run
 *
 * Each process sets up a Laplace problem of its own size and solves it with an iterative and a direct solver
 * for sequential linear algebra. The solution has to satisfy the system and the norm of the solver has to be
 * the 2-norm of the vector of the process.
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/test/laplacian.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>

namespace Dumux::Test {

template<class LinearSolver, class Matrix, class Vector>
void solveOwnSystem(LinearSolver& solver, Matrix A, const Vector& b, const std::string& name)
{
    const auto normDifference = std::abs(solver.norm(b) - b.two_norm());

    // the iterative solvers overwrite the right-hand side
    auto rhs = b;
    Vector x(b.size());
    x = 0.0;
    if (!solver.solve(A, x, rhs))
        DUNE_THROW(Dune::Exception, name << " did not converge");

    auto residual = b;
    A.mmv(x, residual);
    const auto relativeResidual = residual.two_norm()/b.two_norm();
    std::cout << name << " on rank " << Dune::MPIHelper::instance().rank() << ": |b - Ax|/|b| = " << relativeResidual
              << ", |norm(b) - |b|| = " << normDifference << std::endl;

    if (relativeResidual > 1e-8)
        DUNE_THROW(Dune::Exception, name << " does not solve the system of the process");
    if (normDifference > 1e-14*b.two_norm())
        DUNE_THROW(Dune::Exception, "The norm of " << name << " is not the 2-norm of the vector of the process");
}

} // end namespace Dumux::Test

int main(int argc, char* argv[])
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv, "params.input");

    using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    using LATraits = LinearAlgebraTraits<Matrix, Vector>;

    Matrix A;
    setupLaplacian(A, getParam<int>("ProblemSize") + Dune::MPIHelper::instance().rank());
    Vector b(A.N());
    for (std::size_t i = 0; i < b.size(); ++i)
        b[i] = 1.0 + 0.5*i;

    AMGBiCGSTABIstlSolver<SeqLinearSolverTraits, LATraits> iterativeSolver;
    Test::solveOwnSystem(iterativeSolver, A, b, "AMGBiCGSTABIstlSolver");

#if HAVE_UMFPACK
    UMFPackIstlSolver<SeqLinearSolverTraits, LATraits> directSolver;
    Test::solveOwnSystem(directSolver, A, b, "UMFPackIstlSolver");
#endif

    return 0;
}
