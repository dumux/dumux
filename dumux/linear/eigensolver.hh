// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Eigenvalues of a sparse matrix closest to a shift, and their eigenvectors
 */
#ifndef DUMUX_LINEAR_EIGENSOLVER_HH
#define DUMUX_LINEAR_EIGENSOLVER_HH

#if DUMUX_HAVE_SPECTRA && HAVE_UMFPACK

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/istl/foreach.hh>

#include <Eigen/Core>
#include <Spectra/GenEigsRealShiftSolver.h>

#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>

namespace Dumux::Detail::EigenSolver {

/*!
 * \brief The operator \f$ (\mathbf{A} - \sigma \mathbf{I})^{-1} \f$ in the form required by Spectra,
 *        applied with a sparse LU factorization
 */
template<class Matrix, class Vector>
class ShiftInverseOperator
{
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraits<Matrix, Vector>>;

public:
    using Scalar = double;

    ShiftInverseOperator(const Matrix& A, const Vector& blockStructure, Scalar shift)
    : shiftedMatrix_(std::make_shared<Matrix>(A))
    , x_(blockStructure)
    , b_(blockStructure)
    {
        Dune::flatMatrixForEach(*shiftedMatrix_, [&](auto& entry, std::size_t row, std::size_t col)
        {
            if (row == col)
                entry -= shift;
        });
        size_ = Dune::flatVectorForEach(x_, [](auto&&, std::size_t) {});
        solver_.setMatrix(shiftedMatrix_);
    }

    int rows() const { return size_; }
    int cols() const { return size_; }

    // the shift is already part of the factorized matrix
    void set_shift(Scalar) {}

    void perform_op(const Scalar* xIn, Scalar* yOut) const
    {
        Dune::flatVectorForEach(b_, [&](auto& entry, std::size_t i) { entry = xIn[i]; });
        x_ = 0.0;
        solver_.solve(x_, b_);
        Dune::flatVectorForEach(x_, [&](const auto& entry, std::size_t i) { yOut[i] = entry; });
    }

private:
    std::shared_ptr<Matrix> shiftedMatrix_;
    mutable Vector x_, b_;
    mutable LinearSolver solver_;
    int size_;
};

} // end namespace Dumux::Detail::EigenSolver

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Eigenvalues of a sparse matrix closest to a shift \f$ \sigma \f$, and their eigenvectors
 *
 * Arnoldi iteration (Spectra) for the eigenvalues of largest magnitude of
 * \f$ (\mathbf{A} - \sigma \mathbf{I})^{-1} \f$, which belong to the eigenvalues of \f$ \mathbf{A} \f$
 * closest to \f$ \sigma \f$. \f$ \mathbf{A} - \sigma \mathbf{I} \f$ is factorized with UMFPack. The
 * matrix need not be symmetric, so eigenvalues and eigenvectors are complex in general; a real
 * eigenvalue has a real eigenvector. The matrix may be a block or a multi-type block matrix with all
 * diagonal entries stored, and \f$ \sigma \f$ must not be an eigenvalue.
 *
 * Parameters in the group `EigenSolver`:
 * - `Shift` \f$ \sigma \f$ (default 0),
 * - `Tolerance` relative tolerance of the Ritz values (default 1e-10),
 * - `MaxIterations` maximum number of restarts (default 1000),
 * - `NumArnoldiVectors` dimension of the Krylov subspace (default twice the number of requested
 *   eigenvalues plus one and at least 20, at most the size of the matrix).
 */
template<class Matrix, class Vector>
class ShiftInvertEigenSolver
{
    using Scalar = double;

public:
    struct Result
    {
        //! the eigenvalues, ordered by their distance to the shift
        std::vector<std::complex<Scalar>> eigenvalues;
        //! the real and imaginary parts of the eigenvectors
        std::vector<Vector> eigenvectorsReal, eigenvectorsImag;
    };

    explicit ShiftInvertEigenSolver(const std::string& paramGroup = "")
    : shift_(getParamFromGroup<Scalar>(paramGroup, "EigenSolver.Shift", 0.0))
    , tolerance_(getParamFromGroup<Scalar>(paramGroup, "EigenSolver.Tolerance", 1e-10))
    , maxIterations_(getParamFromGroup<int>(paramGroup, "EigenSolver.MaxIterations", 1000))
    , numArnoldiVectors_(getParamFromGroup<int>(paramGroup, "EigenSolver.NumArnoldiVectors", 0))
    {}

    /*!
     * \brief Computes the eigenvalues closest to the shift and their eigenvectors
     * \param A the matrix
     * \param blockStructure a vector with the block structure matching the matrix
     * \param numEigenvalues the number of eigenvalues
     */
    Result solve(const Matrix& A, const Vector& blockStructure, int numEigenvalues) const
    {
        using Operator = Detail::EigenSolver::ShiftInverseOperator<Matrix, Vector>;
        Operator op(A, blockStructure, shift_);

        const int n = op.rows();
        if (numEigenvalues < 1 || numEigenvalues > n - 2)
            DUNE_THROW(Dune::RangeError, "Number of eigenvalues " << numEigenvalues << " not in [1, " << n - 2 << "]");

        const int ncv = std::min(numArnoldiVectors_ > 0 ? numArnoldiVectors_ : std::max(2*numEigenvalues + 1, 20), n);
        Spectra::GenEigsRealShiftSolver<Operator> eigs(op, numEigenvalues, ncv, shift_);
        eigs.init();
        eigs.compute(Spectra::SortRule::LargestMagn, maxIterations_, tolerance_);
        if (eigs.info() != Spectra::CompInfo::Successful)
            DUNE_THROW(NumericalProblem, "Eigenvalue computation did not converge");

        const Eigen::VectorXcd values = eigs.eigenvalues();
        const Eigen::MatrixXcd vectors = eigs.eigenvectors();

        std::vector<int> order(values.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int i, int j)
        { return std::abs(values(i) - shift_) < std::abs(values(j) - shift_); });

        Result result;
        for (int k : order)
        {
            result.eigenvalues.push_back(values(k));
            result.eigenvectorsReal.push_back(blockStructure);
            result.eigenvectorsImag.push_back(blockStructure);
            Dune::flatVectorForEach(result.eigenvectorsReal.back(), [&](auto& entry, std::size_t i) { entry = vectors(i, k).real(); });
            Dune::flatVectorForEach(result.eigenvectorsImag.back(), [&](auto& entry, std::size_t i) { entry = vectors(i, k).imag(); });
        }
        return result;
    }

private:
    Scalar shift_;
    Scalar tolerance_;
    int maxIterations_;
    int numArnoldiVectors_;
};

} // end namespace Dumux

#endif // DUMUX_HAVE_SPECTRA && HAVE_UMFPACK
#endif
