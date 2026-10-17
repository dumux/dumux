// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test of the shift-invert eigenvalue solver against matrices with known eigenvalues
 *
 * - the non-symmetric tridiagonal Toeplitz matrix with diagonal a, subdiagonal b and superdiagonal c
 *   has the eigenvalues \f$ a + 2\sqrt{bc}\cos(k\pi/(n+1)) \f$, \f$ k = 1, \dots, n \f$ (for bc > 0);
 *   it is similar to a symmetric matrix by a diagonal scaling with condition number
 *   \f$ (b/c)^{n/2} \f$, which bounds how accurately the eigenvalues can be computed,
 * - a block-diagonal matrix of 2x2 blocks \f$ ((a_k, -w), (w, a_k)) \f$ has the complex conjugate
 *   eigenvalues \f$ a_k \pm i w \f$,
 * - a multi-type block matrix with two Toeplitz matrices on its diagonal has the eigenvalues of both.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/indices.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/foreach.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/linear/eigensolver.hh>

namespace Dumux {

template<int blockSize>
using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, blockSize, blockSize>>;

template<int blockSize>
using Vector = Dune::BlockVector<Dune::FieldVector<double, blockSize>>;

//! a matrix with the given nonzero pattern per row (column indices), entries zero
template<class M>
M createMatrix(std::size_t rows, std::size_t cols, const std::vector<std::vector<std::size_t>>& pattern)
{
    M A(rows, cols, M::random);
    for (std::size_t i = 0; i < rows; ++i)
        A.setrowsize(i, pattern[i].size());
    A.endrowsizes();
    for (std::size_t i = 0; i < rows; ++i)
        for (auto j : pattern[i])
            A.addindex(i, j);
    A.endindices();
    A = 0.0;
    return A;
}

Matrix<1> toeplitz(std::size_t n, double a, double b, double c)
{
    std::vector<std::vector<std::size_t>> pattern(n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = (i > 0 ? i - 1 : 0); j <= std::min(i + 1, n - 1); ++j)
            pattern[i].push_back(j);

    auto A = createMatrix<Matrix<1>>(n, n, pattern);
    for (std::size_t i = 0; i < n; ++i)
    {
        A[i][i] = a;
        if (i > 0) A[i][i-1] = b;
        if (i + 1 < n) A[i][i+1] = c;
    }
    return A;
}

std::vector<std::complex<double>> toeplitzEigenvalues(std::size_t n, double a, double b, double c)
{
    std::vector<std::complex<double>> values;
    for (std::size_t k = 1; k <= n; ++k)
        values.emplace_back(a + 2.0*std::sqrt(b*c)*std::cos(k*M_PI/(n + 1)));
    return values;
}

//! the given number of exact eigenvalues closest to the shift
std::vector<std::complex<double>> closest(std::vector<std::complex<double>> values, double shift, std::size_t num)
{
    std::sort(values.begin(), values.end(), [&](auto x, auto y) { return std::abs(x - shift) < std::abs(y - shift); });
    values.resize(num);
    return values;
}

// relative residual |A v - lambda v| / |v| of the complex eigenpair
template<class M, class V>
double residual(const M& A, const V& vr, const V& vi, std::complex<double> lambda)
{
    V Avr(vr), Avi(vi);
    A.mv(vr, Avr);
    A.mv(vi, Avi);
    // real part: A vr - (lr vr - li vi), imaginary part: A vi - (lr vi + li vr)
    Avr.axpy(-lambda.real(), vr); Avr.axpy(lambda.imag(), vi);
    Avi.axpy(-lambda.real(), vi); Avi.axpy(-lambda.imag(), vr);
    return std::sqrt(Avr.two_norm2() + Avi.two_norm2())/std::sqrt(vr.two_norm2() + vi.two_norm2());
}

// matches each computed eigenvalue with one of the expected ones (conjugate pairs may come in any order)
void checkEigenvalues(const std::vector<std::complex<double>>& computed,
                      std::vector<std::complex<double>> expected,
                      const std::string& name)
{
    for (const auto& lambda : computed)
    {
        auto it = std::min_element(expected.begin(), expected.end(),
                                   [&](auto x, auto y) { return std::abs(x - lambda) < std::abs(y - lambda); });
        std::cout << name << ": " << lambda << ", deviation " << std::abs(*it - lambda) << std::endl;
        if (std::abs(*it - lambda) > 1e-9)
            DUNE_THROW(Dune::Exception, name << ": eigenvalue " << lambda << " deviates from " << *it);
        expected.erase(it);
    }
}

template<class M, class V, class Result>
void checkResiduals(const M& A, const Result& result, const std::string& name)
{
    for (std::size_t k = 0; k < result.eigenvalues.size(); ++k)
    {
        const auto r = residual(A, result.eigenvectorsReal[k], result.eigenvectorsImag[k], result.eigenvalues[k]);
        if (r > 1e-8)
            DUNE_THROW(Dune::Exception, name << ": eigenpair " << k << " has the residual " << r);
    }
}

void testToeplitz(const std::string& paramGroup, double shift)
{
    const std::size_t n = 100;
    const double a = 2.0, b = -1.05, c = -0.95;
    const auto A = toeplitz(n, a, b, c);

    ShiftInvertEigenSolver<Matrix<1>, Vector<1>> solver(paramGroup);
    const int numEigenvalues = 5;
    const auto result = solver.solve(A, Vector<1>(n), numEigenvalues);

    const std::string name = "Toeplitz, shift " + std::to_string(shift);
    checkEigenvalues(result.eigenvalues, closest(toeplitzEigenvalues(n, a, b, c), shift, numEigenvalues), name);
    checkResiduals<Matrix<1>, Vector<1>>(A, result, name);

    for (std::size_t k = 1; k < result.eigenvalues.size(); ++k)
        if (std::abs(result.eigenvalues[k] - shift) < std::abs(result.eigenvalues[k-1] - shift))
            DUNE_THROW(Dune::Exception, name << ": eigenvalues are not ordered by their distance to the shift");
}

void testComplexPairs()
{
    const std::size_t n = 50;
    const double w = 0.5;
    std::vector<std::vector<std::size_t>> pattern(n);
    for (std::size_t i = 0; i < n; ++i)
        pattern[i] = {i};
    auto A = createMatrix<Matrix<2>>(n, n, pattern);

    std::vector<std::complex<double>> expected;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double ak = 0.1*(i + 1);
        A[i][i] = {{ak, -w}, {w, ak}};
        expected.emplace_back(ak, w);
        expected.emplace_back(ak, -w);
    }

    ShiftInvertEigenSolver<Matrix<2>, Vector<2>> solver;
    const auto result = solver.solve(A, Vector<2>(n), 4);
    checkEigenvalues(result.eigenvalues, closest(expected, 0.0, 4), "complex pairs");
    checkResiduals<Matrix<2>, Vector<2>>(A, result, "complex pairs");
}

void testMultiType()
{
    using namespace Dune::Indices;
    using MTMatrixRow = Dune::MultiTypeBlockVector<Matrix<1>, Matrix<1>>;
    using MTMatrix = Dune::MultiTypeBlockMatrix<MTMatrixRow, MTMatrixRow>;
    using MTVector = Dune::MultiTypeBlockVector<Vector<1>, Vector<1>>;

    const std::size_t n0 = 100, n1 = 80;
    MTMatrix A;
    A[_0][_0] = toeplitz(n0, 2.0, -1.05, -0.95);
    A[_1][_1] = toeplitz(n1, 1.0, -0.45, -0.55);
    A[_0][_1] = createMatrix<Matrix<1>>(n0, n1, std::vector<std::vector<std::size_t>>(n0));
    A[_1][_0] = createMatrix<Matrix<1>>(n1, n0, std::vector<std::vector<std::size_t>>(n1));

    MTVector blockStructure;
    blockStructure[_0].resize(n0);
    blockStructure[_1].resize(n1);

    auto expected = toeplitzEigenvalues(n0, 2.0, -1.05, -0.95);
    const auto second = toeplitzEigenvalues(n1, 1.0, -0.45, -0.55);
    expected.insert(expected.end(), second.begin(), second.end());

    ShiftInvertEigenSolver<MTMatrix, MTVector> solver;
    const auto result = solver.solve(A, blockStructure, 6);
    checkEigenvalues(result.eigenvalues, closest(expected, 0.0, 6), "multi-type");
    checkResiduals<MTMatrix, MTVector>(A, result, "multi-type");
}

} // end namespace Dumux

int main()
{
    Dumux::Parameters::init([](auto& tree) { tree["Shifted.EigenSolver.Shift"] = "1.0"; });
    Dumux::testToeplitz("", 0.0);
    Dumux::testToeplitz("Shifted", 1.0);
    Dumux::testComplexPairs();
    Dumux::testMultiType();
    return 0;
}
