// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Nonlinear
 * \brief The RASPEN preconditioned residual and its exact Jacobian
 */
#ifndef DUMUX_NONLINEAR_PRECONDITIONING_RAS_OPERATOR_HH
#define DUMUX_NONLINEAR_PRECONDITIONING_RAS_OPERATOR_HH

#include <memory>
#include <vector>

#include <dune/istl/bvector.hh>
#include <dune/istl/umfpack.hh>

#include <dumux/assembly/diffmethod.hh>

#include <dumux/nonlinear/preconditioning/localsolver.hh>
#include <dumux/nonlinear/preconditioning/partition.hh>

namespace Dumux::NonlinearPreconditioning {

/*!
 * \ingroup Nonlinear
 * \brief Evaluates \f$ F_{RASPEN}(u) = \sum_i \tilde{R}_i^T C_i(u) \f$ and applies its exact Jacobian
 *
 * The local correction \f$ C_i(u) \f$ is defined by \f$ R_i F(u - R_i^T C_i(u)) = 0 \f$, so
 * \f$ C_i(u) = R_i (u - v_i) \f$ with \f$ v_i \f$ the converged local iterate. Differentiating the
 * defining equation gives \f$ C_i'(u) = A_i^{-1} R_i J(v_i) \f$ with \f$ A_i = R_i J(v_i) R_i^T \f$,
 * evaluated at \f$ v_i \f$ rather than at \f$ u \f$ — this is what makes RASPEN exact where ASPIN is not.
 *
 * Splitting \f$ R_i J(v_i) p = A_i p_i + B_i p_\Gamma \f$ over the subdomain and fringe columns,
 * \f[ J_{RASPEN}(u)\,p = p + \sum_i R_i^T D_i A_i^{-1} B_i p_\Gamma \f]
 * using \f$ \sum_i R_i^T D_i R_i = I \f$. Each application therefore costs one small sparse matrix-vector
 * product and one back-solve per subdomain, with no global matrix ever assembled.
 */
template<class TypeTag, DiffMethod diffMethod = DiffMethod::numeric>
class RestrictedAdditiveSchwarzOperator
{
    using Solver = SubdomainSolver<TypeTag, diffMethod>;

public:
    using SolutionVector = typename Solver::SolutionVector;
    using PrimaryVariables = typename Solver::PrimaryVariables;
    using Problem = typename Solver::Problem;
    using GridGeometry = typename Solver::GridGeometry;
    using GridVariables = typename Solver::GridVariables;
    using LocalVector = typename Solver::LocalVector;
    using SquareMatrix = typename Solver::SquareMatrix;
    using Index = Partition::Index;

    RestrictedAdditiveSchwarzOperator(std::shared_ptr<const Problem> problem,
                   std::shared_ptr<const GridGeometry> gridGeometry,
                   std::shared_ptr<GridVariables> gridVariables,
                   std::shared_ptr<const Partition> partition,
                   const std::vector<std::vector<Index>>& influence)
    : partition_(partition)
    , gridVariables_(gridVariables)
    {
        solvers_.reserve(partition->numSubdomains());
        for (SubdomainIndex i = 0; i < partition->numSubdomains(); ++i)
            solvers_.emplace_back(std::make_unique<Solver>(
                problem, gridGeometry, gridVariables, partition, influence, i));
    }

    using TimeLoop = typename Solver::TimeLoop;

    void setTimeLoop(std::shared_ptr<const TimeLoop> timeLoop)
    {
        for (auto& s : solvers_)
            s->setTimeLoop(timeLoop);
    }

    void setPreviousSolution(const SolutionVector& prevSol)
    {
        for (auto& s : solvers_)
            s->setPreviousSolution(prevSol);
    }

    void setLocalTolerance(double tol)
    {
        for (auto& s : solvers_)
            s->setResidualTolerance(tol);
    }

    void setLocalMaxIterations(std::size_t n)
    {
        for (auto& s : solvers_)
            s->setMaxIterations(n);
    }

    /*!
     * \brief Evaluate the preconditioned residual
     * \param u the current global iterate
     * \param result the preconditioned residual, sized like u
     * \return whether every subdomain solve converged
     *
     * \note Also caches, per subdomain, the Jacobian assembled at the converged local iterate, so a
     *       subsequent applyJacobian() is exact for this u.
     */
    bool evaluate(const SolutionVector& u, SolutionVector& result)
    {
        result = u;
        result = 0.0;

        bool allConverged = true;
        for (SubdomainIndex i = 0; i < partition_->numSubdomains(); ++i)
        {
            SolutionVector v(u);
            const auto localResult = solvers_[i]->solve(v);
            allConverged = allConverged && localResult.converged;

            // the partition of unity is the core indicator, so only the owning subdomain contributes
            for (const auto dof : partition_->core(i))
                result[dof] = u[dof] - v[dof];
        }

        gridVariables_->update(u);
        return allConverged;
    }

    /*!
     * \brief Apply the exact Jacobian of the preconditioned residual
     * \note Requires a preceding evaluate() at the same iterate
     */
    void applyJacobian(const SolutionVector& p, SolutionVector& result) const
    {
        result = p;

        for (SubdomainIndex i = 0; i < partition_->numSubdomains(); ++i)
        {
            const auto& view = solvers_[i]->assembler().jacobian();
            const auto& columns = view.columns();
            const auto& rows = view.rows();

            // B_i p_Gamma: only the columns outside the subdomain contribute
            LocalVector rhs(rows.size());
            rhs = 0.0;
            for (auto row = view.matrix().begin(); row != view.matrix().end(); ++row)
                for (auto col = row->begin(); col != row->end(); ++col)
                {
                    const auto globalCol = columns[col.index()];
                    if (partition_->localIndex(i, globalCol) != Partition::noIndex)
                        continue;
                    col->umv(p[globalCol], rhs[row.index()]);
                }

            LocalVector correction(rows.size());
            correction = 0.0;
            Dune::UMFPack<SquareMatrix> solver(
                const_cast<SquareMatrix&>(solvers_[i]->squareBlock()), 0);
            Dune::InverseOperatorResult solverResult;
            solver.apply(correction, rhs, solverResult);

            for (const auto dof : partition_->core(i))
                result[dof] += correction[partition_->localIndex(i, dof)];
        }
    }

    const Partition& partition() const { return *partition_; }

private:
    std::shared_ptr<const Partition> partition_;
    std::shared_ptr<GridVariables> gridVariables_;
    std::vector<std::unique_ptr<Solver>> solvers_;
};

} // end namespace Dumux::NonlinearPreconditioning

#endif
