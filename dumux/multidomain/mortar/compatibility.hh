// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MortarCoupling
 * \brief Compatibility constraints of the flux mortar on floating subdomains
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_COMPATIBILITY_HH
#define DUMUX_MULTIDOMAIN_MORTAR_COMPATIBILITY_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioner.hh>

#include "projectors.hh"

namespace Dumux::Mortar {

/*!
 * \ingroup MortarCoupling
 * \brief The constraints a flux mortar datum has to satisfy on floating subdomains.
 *
 * A subdomain whose entire boundary carries mortar data is given a flux on all of it, so its
 * problem is a pure Neumann problem: solvable only if the imposed flux balances the sources,
 * and then determined only up to a constant. The constant is an unknown of the coupled
 * problem, since the trace the interface operator reads is a pressure. Written out, the
 * interface problem is
 * \f[ \begin{pmatrix} S & C \\ C^T & 0\end{pmatrix}
 *     \begin{pmatrix} \lambda \\ c\end{pmatrix}
 *   = \begin{pmatrix} b \\ g\end{pmatrix}, \f]
 * where the columns of \f$C\f$ hold the moments \f$\sigma_i\int_{\Gamma_i}\psi_k\f$ of the
 * signed boundary indicator of each floating subdomain, which is both the trace of its
 * constant mode and the functional expressing its solvability.
 *
 * Restricting the iterate to \f$C^T\lambda = g\f$ eliminates \f$c\f$, because the projection
 * annihilates the column space of \f$C\f$. That is the constrained space of \cite Boon2023.
 *
 * \note Only the homogeneous constraint \f$g = 0\f$ is implemented, which is the case of a
 *       floating subdomain without volumetric sources. With sources, \f$g_i\f$ is their
 *       integral over the subdomain and an initial datum satisfying the constraint is needed
 *       in addition to the projection.
 */
template<typename M>
class CompatibilityConstraints
{
    using Scalar = typename M::SolutionVector::field_type;

public:
    using Model = M;
    using SolutionVector = typename M::SolutionVector;

    //! The constraints of the floating subdomains of the given model
    explicit CompatibilityConstraints(std::shared_ptr<Model> model)
    : model_{std::move(model)}
    {
        build_();
    }

    //! Return true if no subdomain floats, in which case nothing has to be constrained
    bool empty() const
    { return columns_.empty(); }

    //! The number of constraints, one per floating subdomain
    std::size_t size() const
    { return columns_.size(); }

    //! Remove from the given vector the part that violates the constraints
    void project(SolutionVector& v) const
    {
        if (columns_.empty())
            return;

        const auto n = columns_.size();
        Dune::DynamicVector<Scalar> rhs(n), beta(n);
        for (std::size_t i = 0; i < n; ++i)
            rhs[i] = columns_[i]*v;
        gramInverse_.mv(rhs, beta);
        for (std::size_t i = 0; i < n; ++i)
            v.axpy(-beta[i], columns_[i]);
    }

    //! The largest constraint violation of the given vector, relative to its size
    Scalar violation(const SolutionVector& v) const
    {
        Scalar result = 0.0;
        const auto norm = v.two_norm();
        for (const auto& column : columns_)
            result = std::max(result, std::abs(column*v)/(column.two_norm()*(norm > 0.0 ? norm : 1.0)));
        return result;
    }

private:
    void build_()
    {
        std::unordered_map<std::size_t, std::vector<Scalar>> integrals;
        model_->decomposition().visitMortars([&] (const auto& mortarPtr) {
            integrals[model_->decomposition().id(*mortarPtr)] = Detail::mortarBasisIntegrals<Scalar>(*mortarPtr);
        });

        std::unordered_map<std::size_t, SolutionVector> columns;
        model_->visitCouplings([&] (const auto& mortar, const auto& solver, std::size_t subDomainId) {
            if (!solver.isFloating())
                return;

            auto& column = columns[subDomainId];
            if (column.size() != model_->numMortarDofs())
            {
                column.resize(model_->numMortarDofs());
                column = 0.0;
            }

            const auto mortarId = model_->decomposition().id(mortar);
            const auto sign = static_cast<Scalar>(model_->orientation(subDomainId, mortarId));
            const auto offset = model_->mortarDofOffset(mortar);
            const auto& mortarIntegrals = integrals.at(mortarId);
            for (std::size_t k = 0; k < mortarIntegrals.size(); ++k)
                column[offset + k] += sign*mortarIntegrals[k];
        });

        if (columns.empty())
            return;

        for (auto& [id, column] : columns)
            columns_.push_back(std::move(column));

        const auto n = columns_.size();
        gramInverse_.resize(n, n);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                gramInverse_[i][j] = columns_[i]*columns_[j];
        gramInverse_.invert();
    }

    std::shared_ptr<Model> model_;
    std::vector<SolutionVector> columns_;
    Dune::DynamicMatrix<Scalar> gramInverse_;
};

/*!
 * \ingroup MortarCoupling
 * \brief Restricts a linear operator to the data satisfying the given constraints.
 */
template<typename Operator, typename Constraints>
class ConstrainedOperator
: public Dune::LinearOperator<typename Operator::domain_type, typename Operator::range_type>
{
    using X = typename Operator::domain_type;
    using Y = typename Operator::range_type;

public:
    //! Restrict the given operator; both arguments must outlive this object
    ConstrainedOperator(Operator& op, const Constraints& constraints)
    : op_{&op}, constraints_{&constraints} {}

    //! Apply the operator to the constrained part of x and constrain the image
    void apply(const X& x, Y& y) const override
    {
        X projected(x);
        constraints_->project(projected);
        op_->apply(projected, y);
        constraints_->project(y);
    }

    //! Add alpha times the image of x to y
    void applyscaleadd(typename X::field_type alpha, const X& x, Y& y) const override
    {
        Y tmp(y);
        apply(x, tmp);
        y.axpy(alpha, tmp);
    }

    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::sequential; }

private:
    Operator* op_;
    const Constraints* constraints_;
};

/*!
 * \ingroup MortarCoupling
 * \brief Restricts a preconditioner to the data satisfying the given constraints.
 */
template<typename X, typename Constraints>
class ConstrainedPreconditioner : public Dune::Preconditioner<X, X>
{
public:
    //! Restrict the given preconditioner; both arguments must outlive this object
    ConstrainedPreconditioner(Dune::Preconditioner<X, X>& prec, const Constraints& constraints)
    : prec_{&prec}, constraints_{&constraints} {}

    void pre(X& x, X& b) override { prec_->pre(x, b); }
    void post(X& x) override { prec_->post(x); }

    //! Apply the preconditioner to the constrained part of the defect and constrain the update
    void apply(X& v, const X& d) override
    {
        X projected(d);
        constraints_->project(projected);
        prec_->apply(v, projected);
        constraints_->project(v);
    }

    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::sequential; }

private:
    Dune::Preconditioner<X, X>* prec_;
    const Constraints* constraints_;
};

} // end namespace Dumux::Mortar

#endif
