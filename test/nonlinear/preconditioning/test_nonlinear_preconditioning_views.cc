// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <cmath>
#include <iostream>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include <dumux/nonlinear/preconditioning/matrixview.hh>
#include <dumux/nonlinear/preconditioning/partition.hh>
#include <dumux/nonlinear/preconditioning/residualview.hh>

namespace {

using Dumux::NonlinearPreconditioning::Partition;
using Dumux::NonlinearPreconditioning::SubdomainIndex;
using Index = Partition::Index;

constexpr int numEq = 2;
using Block = Dune::FieldMatrix<double, numEq, numEq>;
using PrimaryVariables = Dune::FieldVector<double, numEq>;
using GlobalMatrix = Dune::BCRSMatrix<Block>;

int failures = 0;

void check(bool condition, const std::string& what)
{
    if (!condition)
    {
        std::cerr << "FAIL: " << what << std::endl;
        ++failures;
    }
}

std::vector<std::vector<Index>> structuredInfluence(std::size_t nx, std::size_t ny)
{
    std::vector<std::vector<Index>> influence(nx*ny);
    const auto idx = [nx] (std::size_t i, std::size_t j) { return j*nx + i; };
    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i)
        {
            if (i + 1 < nx) { influence[idx(i,j)].push_back(idx(i+1,j)); influence[idx(i+1,j)].push_back(idx(i,j)); }
            if (j + 1 < ny) { influence[idx(i,j)].push_back(idx(i,j+1)); influence[idx(i,j+1)].push_back(idx(i,j)); }
        }
    return influence;
}

//! a deterministic stand-in for the derivative of row's residual with respect to col's primary variables
double entry(Index row, Index col, int i, int j)
{ return 1.0 + std::sin(1.0 + 3.0*double(row) + 7.0*double(col) + 13.0*double(i) + 17.0*double(j)); }

/*!
 * Mimics the cell-centred numeric-differentiation write pattern: visiting element `col` perturbs its own
 * primary variables and writes the resulting derivatives into the rows of `col` itself and of every
 * element whose residual depends on it.
 */
template<class Matrix>
void visitElement(Matrix& A, Index col, const std::vector<std::vector<Index>>& influence)
{
    for (int i = 0; i < numEq; ++i)
        for (int j = 0; j < numEq; ++j)
            A[col][col][i][j] += entry(col, col, i, j);

    for (const auto row : influence[col])
        for (int i = 0; i < numEq; ++i)
            for (int j = 0; j < numEq; ++j)
                A[row][col][i][j] += entry(row, col, i, j);
}

GlobalMatrix assembleGlobal(const std::vector<std::vector<Index>>& influence)
{
    const auto n = influence.size();
    Dune::MatrixIndexSet pattern(n, n);
    for (Index col = 0; col < n; ++col)
    {
        pattern.add(col, col);
        for (const auto row : influence[col])
            pattern.add(row, col);
    }

    GlobalMatrix A;
    pattern.exportIdx(A);
    A = 0.0;

    for (Index col = 0; col < n; ++col)
        visitElement(A, col, influence);

    return A;
}

bool blocksEqual(const Block& a, const Block& b)
{
    for (int i = 0; i < numEq; ++i)
        for (int j = 0; j < numEq; ++j)
            if (std::abs(a[i][j] - b[i][j]) > 1e-14)
                return false;
    return true;
}

/*!
 * The square block: looping only \f$ N_i \f$ must reproduce every entry of the global Jacobian whose row
 * and column both lie in \f$ N_i \f$. This is what the local nonlinear solve needs.
 */
void testSquareBlock(const Partition& p, const std::vector<std::vector<Index>>& influence,
                     const GlobalMatrix& global, SubdomainIndex i)
{
    Dumux::NonlinearPreconditioning::SubdomainMatrixView<Block> view(p, influence, i, false);
    for (const auto col : p.dofs(i))
        visitElement(view, col, influence);

    check(view.numColumnSinkWrites() == 0, "square block: no column-sink writes");

    std::size_t compared = 0;
    for (const auto row : p.dofs(i))
        for (const auto col : p.dofs(i))
        {
            if (!global.exists(row, col))
                continue;
            const auto lr = p.localIndex(i, row), lc = p.localIndex(i, col);
            check(view.matrix().exists(lr, lc), "square block: entry present in view");
            if (view.matrix().exists(lr, lc))
            {
                check(blocksEqual(view.matrix()[lr][lc], global[row][col]),
                      "square block: entry matches the global assembly");
                ++compared;
            }
        }
    check(compared > 0, "square block: comparison was not vacuous");
}

/*!
 * The rectangular block: only the fringe-extended loop populates the columns outside \f$ N_i \f$ that the
 * exact RASPEN Jacobian needs. The same test also confirms the square-only loop leaves them empty, which
 * is why the extra ring is required rather than merely convenient.
 */
void testRectangularBlock(const Partition& p, const std::vector<std::vector<Index>>& influence,
                          const GlobalMatrix& global, SubdomainIndex i)
{
    const auto columns = p.dofsWithFringe(i);
    std::vector<Index> colOf(p.numDofs(), Partition::noIndex);
    for (Index k = 0; k < columns.size(); ++k)
        colOf[columns[k]] = k;

    Dumux::NonlinearPreconditioning::SubdomainMatrixView<Block> view(p, influence, i, true);
    for (const auto col : columns)
        visitElement(view, col, influence);

    check(view.numColumnSinkWrites() == 0, "rectangular block: no column-sink writes");

    std::size_t comparedFringe = 0;
    for (const auto row : p.dofs(i))
        for (const auto col : columns)
        {
            if (!global.exists(row, col))
                continue;
            const auto lr = p.localIndex(i, row), lc = colOf[col];
            check(view.matrix().exists(lr, lc), "rectangular block: entry present in view");
            if (view.matrix().exists(lr, lc))
            {
                check(blocksEqual(view.matrix()[lr][lc], global[row][col]),
                      "rectangular block: entry matches the global assembly");
                if (p.localIndex(i, col) == Partition::noIndex)
                    ++comparedFringe;
            }
        }
    check(comparedFringe > 0, "rectangular block: fringe columns were actually exercised");

    Dumux::NonlinearPreconditioning::SubdomainMatrixView<Block> squareLoop(p, influence, i, true);
    for (const auto col : p.dofs(i))
        visitElement(squareLoop, col, influence);

    bool anyFringeEntryNonZero = false;
    for (const auto row : p.dofs(i))
        for (const auto col : p.fringe(i))
        {
            const auto lr = p.localIndex(i, row), lc = colOf[col];
            if (squareLoop.matrix().exists(lr, lc) && squareLoop.matrix()[lr][lc].infinity_norm() > 1e-14)
                anyFringeEntryNonZero = true;
        }
    check(!anyFringeEntryNonZero,
          "rectangular block: the square-only loop leaves fringe columns empty, so the extra ring is required");
}

void testResidualView(const Partition& p, SubdomainIndex i)
{
    Dumux::NonlinearPreconditioning::SubdomainResidualView<PrimaryVariables> res(p, i);
    check(res.size() == p.dofs(i).size(), "residual view: size matches the subdomain");

    for (Index dof = 0; dof < p.numDofs(); ++dof)
        for (int eq = 0; eq < numEq; ++eq)
            res[dof][eq] = double(dof) + 0.5*double(eq);

    for (const auto dof : p.dofs(i))
        for (int eq = 0; eq < numEq; ++eq)
            check(std::abs(res.vector()[p.localIndex(i, dof)][eq] - (double(dof) + 0.5*double(eq))) < 1e-14,
                  "residual view: in-subdomain entry stored");

    double sum = 0.0;
    for (const auto& block : res.vector())
        sum += block.one_norm();
    check(sum > 0.0, "residual view: comparison was not vacuous");
}

} // end anonymous namespace

int main()
{
    const auto influence = structuredInfluence(6, 5);
    const auto global = assembleGlobal(influence);

    std::vector<SubdomainIndex> cores(30);
    for (std::size_t j = 0; j < 5; ++j)
        for (std::size_t i = 0; i < 6; ++i)
            cores[j*6 + i] = (i < 3 ? 0u : 1u);

    for (std::size_t overlap : {0ul, 1ul, 2ul})
    {
        const Partition p(influence, cores, overlap);
        for (SubdomainIndex i = 0; i < p.numSubdomains(); ++i)
        {
            testSquareBlock(p, influence, global, i);
            testRectangularBlock(p, influence, global, i);
            testResidualView(p, i);
        }
    }

    if (failures > 0)
    {
        std::cerr << failures << " check(s) failed" << std::endl;
        return 1;
    }

    std::cout << "All view checks passed" << std::endl;
    return 0;
}
