//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <iostream>
#include <iomanip>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dune/istl/test/laplacian.hh>
#include <dune/istl/paamg/test/anisotropic.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/linear/linearsolvertraits.hh>

#include <dumux/linear/istlsolverfactorybackend.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>

namespace Dumux::Test {

struct MockGridGeometry
{
    using GridView = Dune::YaspGrid<2>::LeafGridView;
    using DofMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using DiscretizationMethod = DiscretizationMethods::Box;
    static constexpr DiscretizationMethod discMethod{};
};

template<class M, class X, class V>
void solveWithFactory(M& A, X& x, V& b, const std::string& paramGroup)
{
    std::cout << std::endl;

    using LinearSolver = IstlSolverFactoryBackend<LinearSolverTraits<Test::MockGridGeometry>,
                                                  LinearAlgebraTraits<M,V>>;
    LinearSolver solver(paramGroup);

    std::cout << "Solving Laplace problem with " << solver.name() << "\n";
    solver.solve(A, x, b);
    if (!solver.result().converged)
        DUNE_THROW(Dune::Exception, solver.name() << " did not converge!");
}

} // end namespace Dumux::Test

int main(int argc, char* argv[])
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    Parameters::init(argc, argv, "params.input");

    using Vector = Dune::BlockVector<Dune::FieldVector<double, 2>>;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 2, 2>>;

    // create matrix & vectors (we want to solve Ax=b)
    Matrix A; const int numDofs = getParam<int>("ProblemSize");
    setupLaplacian(A, numDofs);

    Vector x(A.N()); Vector b(A.M());
    x = 0; b = 1;

    // AMGBiCGSTABBackend
    {
        std::cout << std::endl;

        const auto testSolverName = "AMGBiCGSTAB";
        using LinearSolver = AMGBiCGSTABIstlSolver<LinearSolverTraits<Test::MockGridGeometry>,
                                                   LinearAlgebraTraits<Matrix,Vector>>;
        LinearSolver solver(testSolverName);

        std::cout << "Solving Laplace problem with " << solver.name() << "\n";
        auto result = solver.solve(A, x, b);
        if (!result.converged)
            DUNE_THROW(Dune::Exception, testSolverName << " did not converge!");

        if (!result)
            DUNE_THROW(Dune::Exception, "Solver result cannot be implicitly converted to bool");
    }

    // IstlSolverFactoryBackend
    Test::solveWithFactory(A, x, b, "AMGCG");
    Test::solveWithFactory(A, x, b, "SSORCG");

    return 0;
}
