// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_EIGENVALUE_PROBLEM_HH
#define DUMUX_EIGENVALUE_PROBLEM_HH

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>

// use Eigen and Spectra to estimate eigenvalues iteratively
// we provide operators that interface Eigen from Dune
#include <Eigen/Core>
#include <Spectra/GenEigsSolver.h>

namespace Dumux::SpectraInterface {

template<class Vector>
void spectraToDuneVector(const double *x_in, Vector& out)
{
    std::size_t startIndex = 0;
    Dune::Hybrid::forEach(out, [&](auto& subVector)
    {
        const auto numEq = std::decay_t<decltype(subVector)>::block_type::size();
        for(std::size_t i = 0; i < subVector.size(); ++i)
            for(std::size_t j = 0; j < numEq; ++j)
                subVector[i][j] = x_in[startIndex + i*numEq + j];
        startIndex += numEq*subVector.size();
    });
}

template<class Vector>
void duneToSpectraVector(const Vector& in, double* y_out)
{
    std::size_t startIndex = 0;
    Dune::Hybrid::forEach(in, [&](const auto& subVector)
    {
        const auto numEq = std::decay_t<decltype(subVector)>::block_type::size();
        for(std::size_t i = 0; i < subVector.size(); ++i)
            for(std::size_t j = 0; j < numEq; ++j)
                y_out[startIndex + i*numEq + j] = subVector[i][j];
        startIndex += numEq*subVector.size();
    });
}

template<class Impl, class Vector>
class DuneMultiTypeOperator
{
public:
    DuneMultiTypeOperator(const Impl& op, const Vector& v)
    : op_(op)
    , v_(v)
    {}

    using Scalar = double;

    int rows() const { return v_.dim(); }
    int cols() const { return v_.dim(); }

    // y_out = OP(x_in)
    void perform_op(const double *x_in, double *y_out) const
    {
        auto x = v_;
        spectraToDuneVector(x_in, x);

        auto out = v_; out = 0.0;
        op_.apply(x, out);

        duneToSpectraVector(out, y_out);
    }
private:
    const Impl& op_;
    const Vector& v_;
};

// Operator: applies K^{-1} * M where M is a diagonal mass matrix given as a vector
// perform_op(x) = K^{-1} * (M * x)
// Used to solve the generalized eigenvalue problem K*v = omega^2 * M*v:
// the largest eigenvalues of K^{-1}*M are 1/omega^2 (lowest frequency modes)
template<class O, class Vector>
class MassWeightedInverseImpl
{
    using M = typename O::matrix_type;
    using MS = decltype(MatrixConverter<M>::multiTypeToBCRSMatrix(std::declval<M>()));
    using VS = decltype(VectorConverter<Vector>::multiTypeToBlockVector(std::declval<Vector>()));
    using LAT = Dumux::LinearAlgebraTraits<M, Vector, MS, VS>;
public:
    MassWeightedInverseImpl(const O& linop, const Vector& massVec)
    : invOp_()
    , massFlat_(massVec.dim())
    {
        matrix_ = linop.getmat();
        invOp_.setMatrix(matrix_);
        duneToSpectraVector(massVec, massFlat_.data());
    }

    // y = K^{-1} * (M * x)
    void apply(const Vector& x, Vector& y) const
    {
        // apply the diagonal mass matrix M: Mx_i = m_i * x_i
        std::vector<double> x_flat(x.dim());
        duneToSpectraVector(x, x_flat.data());
        for (std::size_t i = 0; i < x.dim(); ++i)
            x_flat[i] *= massFlat_[i];

        auto Mx = x;
        spectraToDuneVector(x_flat.data(), Mx);

        // solve K * y = Mx
        invOp_.solve(y, Mx);
    }
private:
    mutable Dumux::UMFPackIstlSolver<SeqLinearSolverTraits, LAT> invOp_;
    M matrix_;
    std::vector<double> massFlat_;
};

} // end namespace Dumux::SpectraInterface

namespace Dumux {

// Solve the generalized eigenvalue problem K*v = omega^2 * M*v.
// Returns the numEigs smallest squared frequencies omega^2 and corresponding eigenvectors.
// massVec is the diagonal of M stored as a solution vector (zero entries for DOFs with no inertia).
// b is used only to provide the correct block structure for the solution vectors.
template<class O, class Vector>
auto smallestNEigenmodes(const O& linop, const Vector& b, const Vector& massVec, int numEigs)
{
    using EigenValues = std::vector<std::complex<typename Vector::field_type>>;
    using EigenVectors = std::vector<Vector>;
    std::tuple<EigenValues, EigenVectors> result;

    // operator: x -> K^{-1} * M * x
    // largest eigenvalues of this operator are 1/omega^2 (lowest frequencies)
    SpectraInterface::MassWeightedInverseImpl<O, Vector> massInvOp(linop, massVec);
    SpectraInterface::DuneMultiTypeOperator<
        SpectraInterface::MassWeightedInverseImpl<O, Vector>, Vector
    > op(massInvOp, b);

    const double tolerance = getParam<double>("EigenSolver.Tolerance", 1e-5);
    const int maxIter = getParam<int>("EigenSolver.MaxIter", 5000);
    const int ncvUser = getParam<int>("EigenSolver.Ncv", 30);
    const int ncvLargest = std::min(std::max(ncvUser, 2*numEigs+2), op.rows());

    Spectra::GenEigsSolver<
        SpectraInterface::DuneMultiTypeOperator<
            SpectraInterface::MassWeightedInverseImpl<O, Vector>, Vector
        >
    > eigs(op, numEigs, ncvLargest);

    eigs.init();
    eigs.compute(Spectra::SortRule::LargestMagn, maxIter, tolerance, Spectra::SortRule::LargestMagn);

    if (eigs.info() != Spectra::CompInfo::Successful)
        DUNE_THROW(Dune::InvalidStateException, "Eigenvalue computation did not converge.");

    const Eigen::VectorXcd levalues = eigs.eigenvalues();
    const Eigen::MatrixXcd ev = eigs.eigenvectors(numEigs);

    auto& [eigenvalues, eigenvectors] = result;
    eigenvalues.resize(numEigs);
    eigenvectors.resize(numEigs);
    for (int k = 0; k < numEigs; ++k)
    {
        // eigenvalue of J^{-1}*M is ±1/omega^2 (sign depends on sign convention of assembled J)
        // return omega^2 = |real(1/eigenvalue)|, always real and positive
        eigenvalues[k] = std::abs(std::real(1.0 / levalues(k)));

        eigenvectors[k] = b;
        std::size_t startIndex = 0;
        Dune::Hybrid::forEach(eigenvectors[k], [&](auto& subVector)
        {
            const auto numEq = std::decay_t<decltype(subVector)>::block_type::size();
            for(std::size_t i = 0; i < subVector.size(); ++i)
                for(std::size_t j = 0; j < numEq; ++j)
                    subVector[i][j] = std::real(ev(startIndex + i*numEq + j, k));
            startIndex += numEq*subVector.size();
        });
    }

    return result;
}

} // end namespace Dumux

#endif
