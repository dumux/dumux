// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test of the two-level p-multigrid preconditioner (PQk -> P1 -> AMG)
 *
 * Assembles reaction-diffusion (symmetric) and convection-diffusion-reaction
 * (non-symmetric) systems with Lagrange elements of order 2 and 3 on the dofs of the
 * FE and hybrid PQk grid geometries and compares Krylov iteration counts for
 * Jacobi, SSOR and p-multigrid preconditioning under grid refinement.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fmatrixev.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/alugrid/grid.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solvers.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/pq2/fegriddiscretization.hh>
#include <dumux/discretization/pq2/fvgridgeometry.hh>
#include <dumux/discretization/pq3/fegriddiscretization.hh>
#include <dumux/discretization/pq3/fvgridgeometry.hh>
#include <dumux/linear/lagrangeprolongation.hh>
#include <dumux/linear/pmultigridpreconditioner.hh>
#include <dumux/linear/pmultigridsmoothers.hh>

namespace Dumux::Test {

template<int bs>
using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, bs, bs>>;
template<int bs>
using Vector = Dune::BlockVector<Dune::FieldVector<double, bs>>;

// Assembles -div(grad u) + beta.grad u + c u = 1 with homogeneous Dirichlet conditions at x_0 = 0.
// The blocks are the scalar entries times the matrix coupling (identity for bs = 1).
template<int bs, class GridGeometry>
void assemble(const GridGeometry& gg, double c, const Dune::FieldVector<double, GridGeometry::GridView::dimensionworld>& beta,
              Matrix<bs>& A, Vector<bs>& b, std::vector<bool>& isDirichlet)
{
    constexpr int dim = GridGeometry::GridView::dimension;
    const auto& gridView = gg.gridView();

    Dune::FieldMatrix<double, bs, bs> coupling(0.0);
    for (int i = 0; i < bs; ++i)
        for (int j = 0; j < bs; ++j)
            coupling[i][j] = (i == j) ? 1.0 : 0.2;

    Dune::MatrixIndexSet pattern(gg.numDofs(), gg.numDofs());
    auto fvGeometry = localView(gg);
    for (const auto& element : elements(gridView))
    {
        fvGeometry.bind(element);
        for (const auto& dofI : localDofs(fvGeometry))
            for (const auto& dofJ : localDofs(fvGeometry))
                pattern.add(dofI.dofIndex(), dofJ.dofIndex());
    }
    pattern.exportIdx(A);
    A = 0.0;
    b.resize(gg.numDofs());
    b = 0.0;
    isDirichlet.assign(gg.numDofs(), false);

    std::vector<Dune::FieldVector<double, 1>> values;
    std::vector<Dune::FieldMatrix<double, 1, dim>> jacobians;
    for (const auto& element : elements(gridView))
    {
        fvGeometry.bind(element);
        const auto geometry = element.geometry();
        const auto& basis = fvGeometry.feLocalBasis();
        const auto order = 2*basis.order();
        const auto& rule = Dune::QuadratureRules<double, dim>::rule(element.type(), order);
        for (const auto& qp : rule)
        {
            basis.evaluateFunction(qp.position(), values);
            basis.evaluateJacobian(qp.position(), jacobians);
            const auto jit = geometry.jacobianInverseTransposed(qp.position());
            std::vector<Dune::FieldVector<double, dim>> grads(values.size());
            for (std::size_t i = 0; i < values.size(); ++i)
                jit.mv(jacobians[i][0], grads[i]);
            const auto w = qp.weight()*geometry.integrationElement(qp.position());

            for (const auto& dofI : localDofs(fvGeometry))
            {
                const auto i = dofI.index();
                for (int k = 0; k < bs; ++k)
                    b[dofI.dofIndex()][k] += w*values[i][0];
                for (const auto& dofJ : localDofs(fvGeometry))
                {
                    const auto j = dofJ.index();
                    const double a = w*(grads[i]*grads[j] + (beta*grads[j])*values[i][0] + c*values[i][0]*values[j][0]);
                    auto block = coupling;
                    block *= a;
                    A[dofI.dofIndex()][dofJ.dofIndex()] += block;
                }
            }
        }

        for (const auto& dof : localDofs(fvGeometry))
            if (ipData(fvGeometry, dof).global()[0] < 1e-10)
                isDirichlet[dof.dofIndex()] = true;
    }

    const bool symmetric = beta.two_norm() == 0.0;
    for (auto row = A.begin(); row != A.end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
            if (isDirichlet[row.index()] || (symmetric && isDirichlet[col.index()]))
            {
                *col = 0.0;
                if (row.index() == col.index())
                    for (int k = 0; k < bs; ++k)
                        (*col)[k][k] = 1.0;
            }
    for (std::size_t i = 0; i < b.size(); ++i)
        if (isDirichlet[i])
            b[i] = 0.0;
}

struct Counts { int jacobi, ssor, pmg; bool symmetric; };

template<int bs, class M, class X>
int solveIterations(M& A, X& x, X b, Dune::Preconditioner<X, X>& prec, bool symmetric, double reduction, int maxIt)
{
    Dune::MatrixAdapter<M, X, X> op(A);
    Dune::SeqScalarProduct<X> sp;
    Dune::InverseOperatorResult res;
    x = 0.0;
    if (symmetric)
    {
        Dune::CGSolver<X> solver(op, sp, prec, reduction, maxIt, 0);
        solver.apply(x, b, res);
    }
    else
    {
        Dune::BiCGSTABSolver<X> solver(op, sp, prec, reduction, maxIt, 0);
        solver.apply(x, b, res);
    }
    return res.converged ? res.iterations : -res.iterations;
}

// checks that prolongation interpolates affine functions and returns the prolongation
template<int bs, class GridGeometry>
auto makeCheckedProlongation(const GridGeometry& gg)
{
    constexpr int dim = GridGeometry::GridView::dimension;
    auto P = makeLagrangeP1Prolongation<bs>(gg);

    Dune::MultipleCodimMultipleGeomTypeMapper<typename GridGeometry::GridView> vertexMapper(gg.gridView(), Dune::mcmgVertexLayout());
    auto affine = [](const auto& x) { double v = 0.3; for (int i = 0; i < dim; ++i) v += (i + 1.0)*x[i]; return v; };

    Vector<bs> xc(P.M()), xf(P.N()), exact(P.N());
    auto fvGeometry = localView(gg);
    for (const auto& element : elements(gg.gridView()))
    {
        fvGeometry.bind(element);
        for (int i = 0; i < element.geometry().corners(); ++i)
            xc[vertexMapper.subIndex(element, i, dim)] = affine(element.geometry().corner(i));
        for (const auto& dof : localDofs(fvGeometry))
            exact[dof.dofIndex()] = affine(ipData(fvGeometry, dof).global());
    }
    P.mv(xc, xf);
    xf -= exact;
    if (xf.infinity_norm() > 1e-12)
        DUNE_THROW(Dune::Exception, "Prolongation does not reproduce affine functions: " << xf.infinity_norm());
    return P;
}

template<int bs, class GridGeometry>
Counts runLevel(const GridGeometry& gg, double c, bool symmetric)
{
    constexpr int dim = GridGeometry::GridView::dimensionworld;
    Dune::FieldVector<double, dim> beta(0.0);
    if (!symmetric)
        for (int i = 0; i < dim; ++i) beta[i] = 5.0*(1.0 - 0.4*i);

    Matrix<bs> A; Vector<bs> b; std::vector<bool> isDirichlet;
    assemble<bs>(gg, c, beta, A, b, isDirichlet);
    auto P = makeCheckedProlongation<bs>(gg);
    Vector<bs> x(b.size());

    const double tol = 1e-8;
    const int maxIt = 3000;

    Counts counts{0, 0, 0, symmetric};
    {
        Dune::SeqJac<Matrix<bs>, Vector<bs>, Vector<bs>> prec(A, 1, 1.0);
        counts.jacobi = solveIterations<bs>(A, x, b, prec, symmetric, tol, maxIt);
    }
    {
        Dune::SeqSSOR<Matrix<bs>, Vector<bs>, Vector<bs>> prec(A, 1, 1.0);
        counts.ssor = solveIterations<bs>(A, x, b, prec, symmetric, tol, maxIt);
    }

    Dune::ParameterTree params;
    params["criterionSymmetric"] = symmetric ? "true" : "false";
    using Factory = PMultigridPreconditionerFactory<typename GridGeometry::GridView, Matrix<bs>>;
    const Factory factory(std::make_shared<const typename Factory::ProlongationMatrix>(std::move(P)));
    const auto op = std::make_shared<Dune::MatrixAdapter<Matrix<bs>, Vector<bs>, Vector<bs>>>(A);
    auto pmg = factory(op, params);
    counts.pmg = solveIterations<bs>(A, x, b, *pmg, symmetric, tol, maxIt);
    if (counts.pmg < 0)
        DUNE_THROW(Dune::Exception, "p-multigrid preconditioned solver did not converge");

    // the residual of the returned solution, independent of the residual norm used inside the solver
    auto residual = b;
    A.mmv(x, residual);
    const auto relativeResidual = residual.two_norm()/b.two_norm();
    if (relativeResidual > 10*tol)
        DUNE_THROW(Dune::Exception, "Relative residual " << relativeResidual << " of the p-multigrid solution");
    return counts;
}

template<class Grid>
inline constexpr bool isBisection = !std::is_same_v<Grid, Dune::YaspGrid<Grid::dimension>>;

template<int bs, class Grid, template<class> class MakeGG>
void runCase(const std::string& name, std::shared_ptr<Grid> grid, int levels, double c, bool symmetric)
{
    std::cout << "\n" << name << ", c = " << c << (symmetric ? ", CG" : ", BiCGSTAB") << ", block size " << bs << "\n"
              << std::setw(10) << "#dofs" << std::setw(10) << "Jacobi" << std::setw(10) << "SSOR" << std::setw(10) << "p-MG\n";
    std::vector<Counts> history;
    for (int level = 0; level < levels; ++level)
    {
        if (level > 0)
            grid->globalRefine(isBisection<Grid> ? Grid::dimension : 1);
        MakeGG<typename Grid::LeafGridView> make;
        auto gg = make(grid->leafGridView());
        const auto counts = runLevel<bs>(*gg, c, symmetric);
        history.push_back(counts);
        std::cout << std::setw(10) << gg->numDofs() << std::setw(10) << counts.jacobi
                  << std::setw(10) << counts.ssor << std::setw(10) << counts.pmg << std::endl;
    }
    const auto& first = history.front();
    const auto& last = history.back();
    if (last.pmg > 2*first.pmg + 2)
        DUNE_THROW(Dune::Exception, "p-multigrid iteration count not bounded: " << first.pmg << " -> " << last.pmg);
    if (!(last.jacobi < 0 || last.pmg < last.jacobi) || !(last.ssor < 0 || last.pmg < last.ssor))
        DUNE_THROW(Dune::Exception, "p-multigrid is not better than Jacobi/SSOR");
}

// checks that shiftSingularValues adds the shift to every singular value and keeps symmetric blocks symmetric
template<int n>
void checkSingularValueShift(const std::string& name, const Dune::FieldMatrix<double, n, n>& D)
{
    const double shift = 0.7;
    auto shifted = D;
    Detail::PMultigrid::shiftSingularValues(shifted, shift);

    auto Q = shifted;
    Q -= D;
    Q /= shift;
    Dune::FieldMatrix<double, n, n> QtQ(0.0), MtM(0.0), StS(0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
            {
                QtQ[i][j] += Q[k][i]*Q[k][j];
                MtM[i][j] += D[k][i]*D[k][j];
                StS[i][j] += shifted[k][i]*shifted[k][j];
            }

    Dune::FieldVector<double, n> before, after;
    Dune::FMatrixHelp::eigenValues(MtM, before);
    Dune::FMatrixHelp::eigenValues(StS, after);
    double orthogonalityError = 0.0, shiftError = 0.0, asymmetry = 0.0;
    bool isSymmetric = true;
    for (int i = 0; i < n; ++i)
    {
        shiftError = std::max(shiftError, std::abs(std::sqrt(std::max(after[i], 0.0)) - std::sqrt(std::max(before[i], 0.0)) - shift));
        for (int j = 0; j < n; ++j)
        {
            orthogonalityError = std::max(orthogonalityError, std::abs(QtQ[i][j] - (i == j ? 1.0 : 0.0)));
            isSymmetric = isSymmetric && D[i][j] == D[j][i];
            asymmetry = std::max(asymmetry, std::abs(shifted[i][j] - shifted[j][i]));
        }
    }

    std::cout << "singular value shift, " << name << ": orthogonality error " << orthogonalityError
              << ", shift error " << shiftError << (isSymmetric ? ", asymmetry " + std::to_string(asymmetry) : "") << std::endl;
    // the singular values of D are computed from D^T D, which resolves small ones only to about the square root of the
    // machine precision
    if (orthogonalityError > 1e-12 || shiftError > 1e-6 || (isSymmetric && asymmetry > 1e-12))
        DUNE_THROW(Dune::Exception, "shiftSingularValues failed for " << name);
}

// checks that the sequential preconditioner of a symmetric matrix is symmetric with each fine smoother
template<int bs, class Grid, template<class> class MakeGG>
void checkSymmetry(const std::string& name, std::shared_ptr<Grid> grid)
{
    MakeGG<typename Grid::LeafGridView> make;
    auto gg = make(grid->leafGridView());
    Matrix<bs> A; Vector<bs> b; std::vector<bool> isDirichlet;
    assemble<bs>(*gg, 1.0, Dune::FieldVector<double, Grid::dimensionworld>(0.0), A, b, isDirichlet);

    using Factory = PMultigridPreconditionerFactory<typename Grid::LeafGridView, Matrix<bs>>;
    const Factory factory(std::make_shared<const typename Factory::ProlongationMatrix>(makeLagrangeP1Prolongation<bs>(*gg)));
    const auto op = std::make_shared<Dune::MatrixAdapter<Matrix<bs>, Vector<bs>, Vector<bs>>>(A);

    Vector<bs> u(A.N()), w(A.N());
    for (std::size_t i = 0; i < A.N(); ++i)
        for (int k = 0; k < bs; ++k)
        {
            u[i][k] = std::sin(1.0 + 0.7*i + 0.3*k);
            w[i][k] = std::cos(0.5 + 1.3*i + 0.1*k);
        }

    const std::vector<std::pair<std::string, int>> smoothers{{"ssor", 0}, {"ilu", 0}, {"ilu", 1}, {"lu", 0}};
    for (const auto& [smoother, iluOrder] : smoothers)
    {
        Dune::ParameterTree params;
        params["pmultigridSmoother"] = smoother;
        params["n"] = std::to_string(iluOrder);
        auto prec = factory(op, params);

        Vector<bs> Mu(A.N()), Mw(A.N()), d(b);
        Mu = 0.0;
        prec->pre(Mu, d);
        prec->apply(Mu, u);
        prec->apply(Mw, w);
        prec->post(Mu);

        const double asymmetry = std::abs(Mu*w - u*Mw)/(Mu.two_norm()*w.two_norm());
        std::cout << "symmetry, " << name << ", block size " << bs << ", smoother " << smoother
                  << (smoother == "ilu" ? "(" + std::to_string(iluOrder) + ")" : "") << ": " << asymmetry << std::endl;
        if (asymmetry > 1e-10)
            DUNE_THROW(Dune::Exception, "p-multigrid preconditioner of a symmetric matrix is not symmetric with smoother " << smoother);
    }
}

template<class GV> struct MakePQ2FE { auto operator()(const GV& gv) const { return std::make_shared<Dumux::Experimental::PQ2FEGridDiscretization<double, GV>>(gv); } };
template<class GV> struct MakePQ3FE { auto operator()(const GV& gv) const { return std::make_shared<Dumux::Experimental::PQ3FEGridDiscretization<double, GV>>(gv); } };
template<class GV> struct MakePQ2Hybrid { auto operator()(const GV& gv) const { return std::make_shared<Dumux::PQ2FVGridGeometry<double, GV, false>>(gv); } };
template<class GV> struct MakePQ3Hybrid { auto operator()(const GV& gv) const { return std::make_shared<Dumux::PQ3FVGridGeometry<double, GV, false>>(gv); } };

template<int dim>
auto makeSimplexGrid(int n)
{
    using Grid = Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>;
    Dune::FieldVector<double, dim> lower(0.0), upper(1.0);
    std::array<unsigned int, dim> cells; cells.fill(n);
    return std::shared_ptr<Grid>(Dune::StructuredGridFactory<Grid>::createSimplexGrid(lower, upper, cells));
}

template<int dim>
auto makeCubeGrid(int n)
{
    Dune::FieldVector<double, dim> upper(1.0);
    std::array<int, dim> cells; cells.fill(n);
    return std::make_shared<Dune::YaspGrid<dim>>(upper, cells);
}

} // end namespace Dumux::Test

int main(int argc, char* argv[])
{
    using namespace Dumux;
    using namespace Dumux::Test;
    Dumux::initialize(argc, argv);

    checkSingularValueShift<1>("zero 1x1", {{0.0}});
    checkSingularValueShift<1>("negative 1x1", {{-2.0}});
    checkSingularValueShift<2>("zero 2x2", {{0.0, 0.0}, {0.0, 0.0}});
    checkSingularValueShift<2>("nonsymmetric singular 2x2", {{0.0, 1.0}, {0.0, 0.0}});
    checkSingularValueShift<2>("symmetric singular 2x2", {{1.0, 0.0}, {0.0, 0.0}});
    checkSingularValueShift<2>("symmetric indefinite 2x2", {{1.0, 0.0}, {0.0, -2.0}});
    checkSingularValueShift<2>("nonsymmetric 2x2", {{0.1, 1.0}, {0.0, 0.2}});
    checkSingularValueShift<3>("nonsymmetric rank 1 3x3", {{0.0, 1.0, 1.0}, {0.0, 2.0, 2.0}, {0.0, 0.0, 0.0}});
    checkSingularValueShift<3>("nonsymmetric rank 2 3x3", {{1.0, 2.0, 3.0}, {0.0, 1.0, 4.0}, {0.0, 0.0, 0.0}});
    checkSingularValueShift<3>("symmetric singular 3x3", {{1.0, -1.0, 0.0}, {-1.0, 2.0, -1.0}, {0.0, -1.0, 1.0}});
    checkSingularValueShift<3>("nonsymmetric 3x3", {{2.0, -1.0, 0.5}, {0.3, -0.1, 1.0}, {-0.7, 0.2, 0.05}});

#if HAVE_SUITESPARSE_UMFPACK
    // the exact LU smoother and a direct solve on the coarsest AMG level need UMFPack
    checkSymmetry<1, Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>, MakePQ2FE>("PQ2 FE, 2D simplex", makeSimplexGrid<2>(8));
    checkSymmetry<2, Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>, MakePQ2FE>("PQ2 FE, 2D simplex", makeSimplexGrid<2>(8));
    checkSymmetry<1, Dune::YaspGrid<2>, MakePQ3FE>("PQ3 FE, 2D cube", makeCubeGrid<2>(8));
#endif

    for (const double c : {1.0, 1000.0})
    {
        runCase<1, Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>, MakePQ2FE>("PQ2 FE, 2D simplex", makeSimplexGrid<2>(4), 4, c, true);
        runCase<1, Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>, MakePQ3FE>("PQ3 FE, 2D simplex", makeSimplexGrid<2>(4), 4, c, true);
        runCase<1, Dune::YaspGrid<2>, MakePQ2FE>("PQ2 FE, 2D cube", makeCubeGrid<2>(4), 4, c, true);
        runCase<1, Dune::YaspGrid<2>, MakePQ3FE>("PQ3 FE, 2D cube", makeCubeGrid<2>(4), 4, c, true);
        runCase<1, Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>, MakePQ2FE>("PQ2 FE, 3D simplex", makeSimplexGrid<3>(2), 3, c, true);
        runCase<1, Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>, MakePQ3FE>("PQ3 FE, 3D simplex", makeSimplexGrid<3>(2), 3, c, true);
        runCase<1, Dune::YaspGrid<3>, MakePQ2FE>("PQ2 FE, 3D cube", makeCubeGrid<3>(2), 3, c, true);
    }

    runCase<2, Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>, MakePQ2FE>("PQ2 FE, 2D simplex", makeSimplexGrid<2>(4), 3, 1.0, true);

    runCase<1, Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>, MakePQ2Hybrid>("PQ2 hybrid, 2D simplex", makeSimplexGrid<2>(4), 4, 1.0, false);
    runCase<1, Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>, MakePQ2Hybrid>("PQ2 hybrid, 3D simplex", makeSimplexGrid<3>(2), 3, 1.0, false);
    runCase<1, Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>, MakePQ3Hybrid>("PQ3 hybrid, 2D simplex", makeSimplexGrid<2>(4), 4, 1.0, false);

    std::cout << "\nAll p-multigrid tests passed." << std::endl;
    return 0;
}
