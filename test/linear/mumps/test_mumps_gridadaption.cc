// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Test for the MUMPS solver backend on a distributed grid: a cell-centered two-point flux
 *        system is solved, the grid is refined, the solver is updated with
 *        updateAfterGridAdaption and the refined system is solved, for one domain and for two
 *        coupled subdomains. Also checks the constructor without grid view (sequential runs only)
 *        and the constructor with grid view for a grid that cannot communicate.
 *
 * The two-point flux approximation on a Cartesian grid reproduces linear solutions exactly, so
 * the discrete solution equals the exact one up to round-off.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
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

template<class GlobalPosition>
double exactSolution(const GlobalPosition& x, int domainIdx)
{
    const double y = x.size() > 1 ? x[1] : 0.0;
    return domainIdx == 0 ? 1.0 + x[0] + 2.0*y : 3.0 - 2.0*x[0] + 0.5*y;
}

/*!
 * \brief Two-point flux discretization of -div(grad u_d) + c (u_d - u_other) = c (u_d - u_other)|exact
 *        with Dirichlet conditions from the exact solution. Returns the diagonal block and the
 *        right-hand side of subdomain d; the coupling block is -c|V| on the diagonal.
 */
template<class GridGeometry>
void assemble(const GridGeometry& gridGeometry, int domainIdx, double c,
              Matrix& A, Matrix& coupling, Vector& b)
{
    const auto& gridView = gridGeometry.gridView();
    const auto& mapper = gridGeometry.elementMapper();
    const std::size_t numDofs = gridView.size(0);

    Dune::MatrixIndexSet pattern(numDofs, numDofs);
    Dune::MatrixIndexSet couplingPattern(numDofs, numDofs);
    for (const auto& element : elements(gridView))
    {
        const auto i = mapper.index(element);
        pattern.add(i, i);
        couplingPattern.add(i, i);
        for (const auto& intersection : intersections(gridView, element))
            if (intersection.neighbor())
                pattern.add(i, mapper.index(intersection.outside()));
    }
    pattern.exportIdx(A);
    couplingPattern.exportIdx(coupling);
    A = 0.0;
    coupling = 0.0;
    b.resize(numDofs);
    b = 0.0;

    for (const auto& element : elements(gridView))
    {
        const auto i = mapper.index(element);
        const auto geometry = element.geometry();
        const auto center = geometry.center();
        for (const auto& intersection : intersections(gridView, element))
        {
            const auto faceGeometry = intersection.geometry();
            if (intersection.neighbor())
            {
                const auto j = mapper.index(intersection.outside());
                const double t = faceGeometry.volume()/(intersection.outside().geometry().center() - center).two_norm();
                A[i][i] += t;
                A[i][j] -= t;
            }
            else if (intersection.boundary())
            {
                const double t = faceGeometry.volume()/(faceGeometry.center() - center).two_norm();
                A[i][i] += t;
                b[i] += t*exactSolution(faceGeometry.center(), domainIdx);
            }
        }

        const double reaction = c*geometry.volume();
        A[i][i] += reaction;
        coupling[i][i] = -reaction;
        b[i] += reaction*(exactSolution(center, domainIdx) - exactSolution(center, 1 - domainIdx));
    }
}

template<class GridGeometry>
void checkSolution(const GridGeometry& gridGeometry, const Vector& x, int domainIdx, const std::string& name)
{
    double error = 0.0;
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        const auto i = gridGeometry.elementMapper().index(element);
        error = std::max(error, std::abs(x[i] - exactSolution(element.geometry().center(), domainIdx)));
    }
    error = gridGeometry.gridView().comm().max(error);

    if (gridGeometry.gridView().comm().rank() == 0)
        std::cout << name << ": max error " << error << std::endl;
    if (!(error < 1e-11))
        DUNE_THROW(Dune::Exception, name << ": wrong solution, max error " << error);
}

template<class Solver, class GridGeometry>
void solveSingleDomain(Solver& solver, const GridGeometry& gridGeometry, const std::string& name)
{
    Matrix A, coupling;
    Vector b;
    assemble(gridGeometry, 0, 0.0, A, coupling, b);
    Vector x(b.size());
    x = 0.0;
    solver.solve(A, x, b);
    if (!solver.result().converged)
        DUNE_THROW(Dune::Exception, name << ": the solver did not report convergence");
    checkSolution(gridGeometry, x, 0, name);
}

template<class Solver, class GridGeometry>
void solveTwoDomains(Solver& solver, const GridGeometry& gridGeometry, const std::string& name)
{
    using namespace Dune::Indices;
    MTMatrix A;
    MTVector x, b;
    assemble(gridGeometry, 0, 10.0, A[_0][_0], A[_0][_1], b[_0]);
    assemble(gridGeometry, 1, 10.0, A[_1][_1], A[_1][_0], b[_1]);
    x = b;
    x = 0.0;
    solver.solve(A, x, b);
    checkSolution(gridGeometry, x[_0], 0, name + " (subdomain 0)");
    checkSolution(gridGeometry, x[_1], 1, name + " (subdomain 1)");
}

} // end namespace Dumux::Test

int main(int argc, char* argv[])
{
    using namespace Dumux;
    using namespace Dumux::Test;
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    using Grid = Dune::YaspGrid<2>;
    Grid grid({1.0, 1.0}, {8, 8});
    using GridGeometry = CCTpfaFVGridGeometry<Grid::LeafGridView>;
    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());

    using SolverTraits = LinearSolverTraits<GridGeometry>;
    DirectSolverMumps<SolverTraits, LinearAlgebraTraits> solver(gridGeometry->gridView(), gridGeometry->elementMapper());
    DirectSolverMumps<SolverTraits, MTLinearAlgebraTraits> mtSolver(std::make_tuple(gridGeometry, gridGeometry));
    solveSingleDomain(solver, *gridGeometry, "coarse grid");
    solveTwoDomains(mtSolver, *gridGeometry, "coarse grid, two subdomains");

    const std::size_t numCoarseDofs = gridGeometry->gridView().size(0);
    grid.globalRefine(1);
    gridGeometry->update(grid.leafGridView());

    // on a distributed grid, the solver needs the new dof distribution
    if (mpiHelper.size() > 1)
    {
        bool detected = false;
        try { solveSingleDomain(solver, *gridGeometry, "refined grid without update"); }
        catch (const Dune::InvalidStateException&) { detected = true; }
        if (!detected)
            DUNE_THROW(Dune::Exception, "Solving after grid adaption without updateAfterGridAdaption did not throw");

        // the failed setMatrix also dropped the factorization of the coarse system
        Vector x(numCoarseDofs), b(numCoarseDofs);
        b = 1.0;
        detected = false;
        try { solver.solve(x, b); }
        catch (const Dune::InvalidStateException&) { detected = true; }
        if (!detected)
            DUNE_THROW(Dune::Exception, "Solving with the factorization of the coarse system after a failed setMatrix did not throw");
    }

    solver.updateAfterGridAdaption(gridGeometry->gridView(), gridGeometry->elementMapper());
    mtSolver.updateAfterGridAdaption(std::make_tuple(gridGeometry, gridGeometry));
    solveSingleDomain(solver, *gridGeometry, "refined grid");
    solveTwoDomains(mtSolver, *gridGeometry, "refined grid, two subdomains");

    using SequentialSolver = DirectSolverMumps<SolverTraits, LinearAlgebraTraits>;
    if (mpiHelper.size() == 1)
    {
        SequentialSolver sequentialSolver;
        solveSingleDomain(sequentialSolver, *gridGeometry, "constructor without grid view");
    }
    else
    {
        bool detected = false;
        try { SequentialSolver sequentialSolver; }
        catch (const Dune::InvalidStateException&) { detected = true; }
        if (!detected)
            DUNE_THROW(Dune::Exception, "Constructing the solver without grid view in a parallel run did not throw");
    }

    // every process solves its own system on a grid that cannot communicate
    Dune::OneDGrid lineGrid(16, 0.0, 1.0);
    using LineGridGeometry = CCTpfaFVGridGeometry<Dune::OneDGrid::LeafGridView>;
    LineGridGeometry lineGridGeometry(lineGrid.leafGridView());
    DirectSolverMumps<LinearSolverTraits<LineGridGeometry>, LinearAlgebraTraits>
        lineSolver(lineGridGeometry.gridView(), lineGridGeometry.elementMapper());
    solveSingleDomain(lineSolver, lineGridGeometry, "grid that cannot communicate");

    return 0;
}
