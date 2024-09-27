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
#include <dune/common/classname.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>

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

    std::cout << "Solving Laplace problem with " << solver.name() << " (matrix type: " << Dune::className<M>() << ")\n";
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
    Matrix A;
    const int numDofs = getParam<int>("ProblemSize");
    setupLaplacian(A, numDofs);

    Vector x(A.N()); Vector b(A.M());

    // AMGBiCGSTABBackend
    {
        x = 0; b = 1;
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
    {
        x = 0; b = 1;
        Test::solveWithFactory(A, x, b, "AMGCG");

        x = 0; b = 1;
        Test::solveWithFactory(A, x, b, "SSORCG");
    }

    // MultiTypeBlockMatrix and IstlSolverFactoryBackend
    {
        using MTVector = Dune::MultiTypeBlockVector<Vector, Vector, Vector>;
        using MTMatrixRow = Dune::MultiTypeBlockVector<Matrix, Matrix, Matrix>;
        using MTMatrix = Dune::MultiTypeBlockMatrix<MTMatrixRow, MTMatrixRow, MTMatrixRow>;
        MTMatrix MTA; MTVector MTx, MTb;

        // reset matrix & vectors
        setupLaplacian(A, numDofs);
        x = 0; b = 1;

        Dune::Hybrid::forEach(std::make_index_sequence<MTA.N()>{}, [&](auto i){
            MTA[i][i] = A;
            MTx[i] = x;
            MTb[i] = b;

            Dune::Hybrid::forEach(std::make_index_sequence<MTA.M()>{}, [&](auto j){
                if constexpr (i != j)
                {
                    MTA[i][j] = A; // set some pattern
                    MTA[i][j] = 0; // set all entries to zero to solve uncoupled problems
                }
            });
        });

        Test::solveWithFactory(MTA, MTx, MTb, "SSORBiCGSTAB");
    }

    return 0;
}
