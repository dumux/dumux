// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Test for the MUMPS solver backend with dofs on the element faces, on a non-overlapping
 *        and an overlapping decomposition, for one domain, for two coupled subdomains, for the
 *        two coupled fields as one system with two unknowns per dof, and for two subdomains on
 *        two grids with different resolutions and decompositions
 *
 * Every element adds its contributions to its own faces only. In the non-overlapping
 * decomposition, the rows of the faces on the process border are therefore split among the
 * processes sharing the face. The discretization reproduces linear solutions exactly, so the
 * discrete solution equals the exact one up to round-off.
 */
#include <config.h>

#include <algorithm>
#include <bitset>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dumux/discretization/facecentered/diamond/fvgridgeometry.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/mumpssolver.hh>

namespace Dumux::Test {

using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
using MTMatrix = Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<Matrix, Matrix>,
                                            Dune::MultiTypeBlockVector<Matrix, Matrix>>;
using MTVector = Dune::MultiTypeBlockVector<Vector, Vector>;

using BlockMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 2, 2>>;
using BlockVector = Dune::BlockVector<Dune::FieldVector<double, 2>>;

struct LinearAlgebraTraits { using Matrix = Test::Matrix; using Vector = Test::Vector; };
struct MTLinearAlgebraTraits { using Matrix = MTMatrix; using Vector = MTVector; };
struct BlockLinearAlgebraTraits { using Matrix = BlockMatrix; using Vector = BlockVector; };

template<class GlobalPosition>
double exactSolution(const GlobalPosition& x, int domainIdx)
{ return domainIdx == 0 ? 1.0 + x[0] + 2.0*x[1] : 3.0 - 2.0*x[0] + 0.5*x[1]; }

/*!
 * \brief Discretization of -div(grad u_d) + c (u_d - u_other) = c (u_d - u_other)|exact with
 *        Dirichlet conditions from the exact solution on the boundary faces. Returns the diagonal
 *        block and the right-hand side of subdomain d; the coupling block is -c|V|/n_faces on the
 *        diagonal.
 *
 * The gradient on an element is reconstructed from the face values by the divergence theorem,
 * grad u = sum_j |F_j| n_j u_j / |V|, which is exact for linear u. The flux of face i of an element
 * is |F_i| n_i * grad u; the fluxes of the two elements sharing a face cancel for linear u.
 */
template<class GridGeometry>
void assemble(const GridGeometry& gridGeometry, int domainIdx, double c,
              Matrix& A, Matrix& coupling, Vector& b)
{
    const auto& gridView = gridGeometry.gridView();
    const auto& mapper = gridGeometry.dofMapper();
    const std::size_t numDofs = mapper.size();
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Geometry::GlobalCoordinate;

    Dune::MatrixIndexSet pattern(numDofs, numDofs);
    Dune::MatrixIndexSet couplingPattern(numDofs, numDofs);
    std::vector<bool> isBoundary(numDofs, false);
    std::vector<GlobalPosition> dofPosition(numDofs);
    for (const auto& element : elements(gridView))
    {
        for (const auto& intersection : intersections(gridView, element))
        {
            const auto i = mapper.subIndex(element, intersection.indexInInside(), 1);
            dofPosition[i] = intersection.geometry().center();
            if (intersection.boundary())
                isBoundary[i] = true;
            couplingPattern.add(i, i);
            for (const auto& other : intersections(gridView, element))
                pattern.add(i, mapper.subIndex(element, other.indexInInside(), 1));
        }
    }
    pattern.exportIdx(A);
    couplingPattern.exportIdx(coupling);
    A = 0.0;
    coupling = 0.0;
    b.resize(numDofs);
    b = 0.0;

    for (const auto& element : elements(gridView))
    {
        const double volume = element.geometry().volume();
        std::vector<std::size_t> dofs;
        std::vector<GlobalPosition> gradientWeights;
        for (const auto& intersection : intersections(gridView, element))
        {
            dofs.push_back(mapper.subIndex(element, intersection.indexInInside(), 1));
            auto weight = intersection.centerUnitOuterNormal();
            weight *= intersection.geometry().volume()/volume;
            gradientWeights.push_back(weight);
        }

        const double reaction = c*volume/dofs.size();
        for (std::size_t localI = 0; localI < dofs.size(); ++localI)
        {
            const auto i = dofs[localI];
            if (isBoundary[i])
                continue;

            for (std::size_t localJ = 0; localJ < dofs.size(); ++localJ)
                A[i][dofs[localJ]] += volume*(gradientWeights[localI]*gradientWeights[localJ]);

            A[i][i] += reaction;
            coupling[i][i] -= reaction;
            b[i] += reaction*(exactSolution(dofPosition[i], domainIdx) - exactSolution(dofPosition[i], 1 - domainIdx));
        }
    }

    for (std::size_t i = 0; i < numDofs; ++i)
    {
        if (isBoundary[i])
        {
            A[i][i] = 1.0;
            b[i] = exactSolution(dofPosition[i], domainIdx);
        }
    }
}

template<class GridGeometry>
void checkSolution(const GridGeometry& gridGeometry, const Vector& x, int domainIdx, const std::string& name)
{
    const auto& gridView = gridGeometry.gridView();
    double error = 0.0;
    for (const auto& element : elements(gridView))
    {
        for (const auto& intersection : intersections(gridView, element))
        {
            const auto i = gridGeometry.dofMapper().subIndex(element, intersection.indexInInside(), 1);
            error = std::max(error, std::abs(x[i] - exactSolution(intersection.geometry().center(), domainIdx)));
        }
    }
    error = gridView.comm().max(error);

    if (gridView.comm().rank() == 0)
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

//! The system of the two coupled subdomains, assembled with one block of two unknowns per dof
template<class Solver, class GridGeometry>
void solveBlockSystem(Solver& solver, const GridGeometry& gridGeometry, const std::string& name)
{
    Matrix A0, A1, coupling0, coupling1;
    Vector b0, b1;
    assemble(gridGeometry, 0, 10.0, A0, coupling0, b0);
    assemble(gridGeometry, 1, 10.0, A1, coupling1, b1);

    BlockMatrix A;
    Dune::MatrixIndexSet pattern(A0.N(), A0.M());
    pattern.import(A0);
    pattern.exportIdx(A);
    A = 0.0;
    for (auto rowIt = A0.begin(); rowIt != A0.end(); ++rowIt)
    {
        for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
        {
            const auto i = rowIt.index();
            const auto j = colIt.index();
            A[i][j][0][0] = (*colIt)[0][0];
            A[i][j][1][1] = A1[i][j][0][0];
        }
    }

    BlockVector b(b0.size());
    for (std::size_t i = 0; i < b.size(); ++i)
    {
        A[i][i][0][1] = coupling0[i][i][0][0];
        A[i][i][1][0] = coupling1[i][i][0][0];
        b[i][0] = b0[i][0];
        b[i][1] = b1[i][0];
    }

    BlockVector x(b.size());
    x = 0.0;
    solver.solve(A, x, b);

    Vector x0(x.size()), x1(x.size());
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        x0[i] = x[i][0];
        x1[i] = x[i][1];
    }
    checkSolution(gridGeometry, x0, 0, name + " (component 0)");
    checkSolution(gridGeometry, x1, 1, name + " (component 1)");
}

void runTest(int overlap)
{
    using Grid = Dune::YaspGrid<2>;
    Grid grid({1.0, 1.0}, {10, 10}, std::bitset<2>(), overlap);
    using GridGeometry = FaceCenteredDiamondFVGridGeometry<Grid::LeafGridView>;
    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());

    using SolverTraits = LinearSolverTraits<GridGeometry>;
    DirectSolverMumps<SolverTraits, LinearAlgebraTraits> solver(gridGeometry->gridView(), gridGeometry->dofMapper());
    DirectSolverMumps<SolverTraits, MTLinearAlgebraTraits> mtSolver(std::make_tuple(gridGeometry, gridGeometry));
    DirectSolverMumps<SolverTraits, BlockLinearAlgebraTraits> blockSolver(gridGeometry->gridView(), gridGeometry->dofMapper());

    const std::string name = "overlap " + std::to_string(overlap);
    solveSingleDomain(solver, *gridGeometry, name);
    solveTwoDomains(mtSolver, *gridGeometry, name + ", two subdomains");
    solveBlockSystem(blockSolver, *gridGeometry, name + ", two unknowns per dof");
}

//! Two subdomains without coupling on two grids, one non-overlapping and one overlapping
void runTwoGridsTest()
{
    using Grid = Dune::YaspGrid<2>;
    Grid grid0({1.0, 1.0}, {10, 10}, std::bitset<2>(), 0);
    Grid grid1({1.0, 1.0}, {7, 13}, std::bitset<2>(), 1);
    using GridGeometry = FaceCenteredDiamondFVGridGeometry<Grid::LeafGridView>;
    auto gridGeometry0 = std::make_shared<GridGeometry>(grid0.leafGridView());
    auto gridGeometry1 = std::make_shared<GridGeometry>(grid1.leafGridView());

    using namespace Dune::Indices;
    MTMatrix A;
    MTVector x, b;
    Matrix sameGridCoupling;
    assemble(*gridGeometry0, 0, 0.0, A[_0][_0], sameGridCoupling, b[_0]);
    assemble(*gridGeometry1, 1, 0.0, A[_1][_1], sameGridCoupling, b[_1]);
    Dune::MatrixIndexSet(b[_0].size(), b[_1].size()).exportIdx(A[_0][_1]);
    Dune::MatrixIndexSet(b[_1].size(), b[_0].size()).exportIdx(A[_1][_0]);
    x = b;
    x = 0.0;

    using SolverTraits = LinearSolverTraits<GridGeometry>;
    DirectSolverMumps<SolverTraits, MTLinearAlgebraTraits> solver(std::make_tuple(gridGeometry0, gridGeometry1));
    solver.solve(A, x, b);
    checkSolution(*gridGeometry0, x[_0], 0, "two grids (subdomain 0)");
    checkSolution(*gridGeometry1, x[_1], 1, "two grids (subdomain 1)");
}

} // end namespace Dumux::Test

int main(int argc, char* argv[])
{
    Dune::MPIHelper::instance(argc, argv);

    Dumux::Test::runTest(0);
    Dumux::Test::runTest(1);
    Dumux::Test::runTwoGridsTest();

    return 0;
}
