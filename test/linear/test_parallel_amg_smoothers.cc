//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <iostream>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/linear/helmholtzoperator.hh>
#include <dumux/linear/preconditioners.hh>

template<class Smoother, class LOP, class Vector>
void testSmoother(LOP& lop, const Vector& b)
{
    using Comm = Dune::Amg::SequentialInformation;
    using Amg = Dune::Amg::AMG<LOP, Vector, Smoother, Comm>;
    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
    using Criterion = Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<typename LOP::matrix_type, Dune::Amg::FirstDiagonal>>;

    Dune::Amg::Parameters params(15, 2000, 1.2, 1.6, Dune::Amg::atOnceAccu);
    Criterion criterion(params);
    SmootherArgs smootherArgs;
    smootherArgs.iterations = 2;
    smootherArgs.relaxationFactor = 0.8;
    auto comm = std::make_shared<Comm>();
    auto amg = std::make_shared<Amg>(lop, criterion, smootherArgs, *comm);
    Dune::CGSolver<Vector> solver(lop, *amg, 1e-15, 200, 1);
    Dune::InverseOperatorResult result;

    Dune::Timer timer;

    auto x = b; x = 0.0; auto bTmp = b;
    solver.apply(x, bTmp, result);
    if (!result.converged)
        DUNE_THROW(Dune::Exception, "Solver did not converge!");

    std::cout << "Solver took " << timer.elapsed() << " seconds" << std::endl;

    using std::abs;
    if (abs(x.two_norm2() - x.size()) > 1e-10*x.size())
        DUNE_THROW(Dune::Exception,
            "Wrong result: " << x.two_norm2()
            << " (expected " << x.size() << ")"
            << " difference " << abs(x.two_norm2() - x.size())
            << " tolerance " << 1e-10*x.size()
        );
}

int main(int argc, char* argv[])
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);

    Parameters::init([](Dune::ParameterTree& p){
        p["Problem.Name"] = "problem";
    });

    using Grid = Dune::YaspGrid<2>;
    using Pos = Dune::FieldVector<typename Grid::ctype, 2>;
    std::array<int, 2> cells{{300, 300}};
    auto grid = std::make_shared<Grid>(Pos(1.0), cells);

    Dune::Timer timer;
    auto helmholtzMatrix = makeHelmholtzMatrix<
        Dumux::Properties::TTag::CCTpfaModel
    >(grid->leafGridView(), 1.0, 1.0);
    std::cout << "Operator assembly took " << timer.elapsed() << " seconds" << std::endl;

    using Matrix = std::decay_t<decltype(helmholtzMatrix)>::element_type;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    Vector x(helmholtzMatrix->N());
    Vector b(helmholtzMatrix->M());
    x = 1.0;

    // test the other interface too
    timer.reset();
    auto lop = makeHelmholtzLinearOperator<
        Dune::MatrixAdapter<Matrix, Vector, Vector>,
        Dumux::Properties::TTag::CCTpfaModel
    >(grid->leafGridView(), 1.0, 1.0);
    lop->apply(x, b);
    std::cout << "Operator assembly took " << timer.elapsed() << " seconds" << std::endl;

    testSmoother<Dumux::ParMTSSOR<Matrix, Vector, Vector>>(*lop, b);
    testSmoother<Dumux::ParMTSOR<Matrix, Vector, Vector>>(*lop, b);
    testSmoother<Dumux::ParMTJac<Matrix, Vector, Vector>>(*lop, b);

    return 0;
}
