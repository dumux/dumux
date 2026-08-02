// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Nonlinear
 * \brief The local nonlinear solve defining a subdomain's correction
 */
#ifndef DUMUX_NONLINEAR_PRECONDITIONING_LOCAL_SOLVER_HH
#define DUMUX_NONLINEAR_PRECONDITIONING_LOCAL_SOLVER_HH

#include <cmath>
#include <functional>
#include <memory>
#include <vector>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/umfpack.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/linear/dunevectors.hh>

#include <dumux/nonlinear/preconditioning/assembler.hh>
#include <dumux/nonlinear/preconditioning/localupdate.hh>
#include <dumux/nonlinear/preconditioning/pvswitch.hh>
#include <dumux/nonlinear/preconditioning/partition.hh>

namespace Dumux::NonlinearPreconditioning {

struct LocalSolveResult
{
    bool converged = false;
    std::size_t iterations = 0;
    std::size_t switched = 0;   //!< degrees of freedom that changed phase state during this solve
    double residualNorm = 0.0;
};

/*!
 * \ingroup Nonlinear
 * \brief Solves \f$ R_i F(u - R_i^T C_i(u)) = 0 \f$ for one subdomain, exterior frozen
 *
 * The iteration deliberately assembles *before* testing convergence, and tests the residual rather than
 * the update. A conventional Newton loop tests the shift after the update and therefore exits holding a
 * Jacobian assembled one step before the iterate it reports; reusing that as \f$ J(v_i) \f$ silently
 * degrades RASPEN's exact Jacobian to an approximation of the same character as ASPIN's. Here the matrix
 * left behind by a converged solve is assembled at exactly the returned iterate.
 *
 * The cost of that guarantee is one extra assembly per solve relative to a shift-criterion loop.
 *
 * The assembler carries the fringe columns throughout, so a converged solve leaves behind both the square
 * block \f$ A_i \f$ used for the local linear solves and the rectangular \f$ R_i J(v_i) \f$ the outer
 * RASPEN operator needs.
 */
template<class TypeTag, DiffMethod diffMethod = DiffMethod::numeric>
class SubdomainSolver
{
    using Assembler = SubdomainAssembler<TypeTag, diffMethod>;
    using MatrixBlock = typename Assembler::JacobianMatrix::MatrixBlock;

public:
    using SolutionVector = typename Assembler::SolutionVector;
    using PrimaryVariables = typename Assembler::PrimaryVariables;
    using Problem = typename Assembler::Problem;
    using GridGeometry = typename Assembler::GridGeometry;
    using GridVariables = typename Assembler::GridVariables;
    using SquareMatrix = Dune::BCRSMatrix<MatrixBlock>;
    // linear algebra needs the plain vector type: with a primary variable switch the solution
    // vector's blocks additionally carry a phase state, which ISTL solvers cannot consume
    using LocalVector = typename Dumux::Detail::NativeDuneVectorType<Dune::BlockVector<PrimaryVariables>>::type;
    using Index = Partition::Index;

    SubdomainSolver(std::shared_ptr<const Problem> problem,
                    std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<GridVariables> gridVariables,
                    std::shared_ptr<const Partition> partition,
                    const std::vector<std::vector<Index>>& influence,
                    SubdomainIndex subdomain)
    : assembler_(problem, gridGeometry, gridVariables, partition, influence, subdomain, true)
    , gridVariables_(gridVariables)
    , gridGeometry_(gridGeometry)
    , partition_(partition)
    , subdomain_(subdomain)
    , updateRange_(partition->dofsWithFringe(subdomain))
    {
        const auto& view = assembler_.jacobian();
        Dune::MatrixIndexSet pattern(view.N(), view.N());
        for (auto row = view.matrix().begin(); row != view.matrix().end(); ++row)
            for (auto col = row->begin(); col != row->end(); ++col)
                if (partition->localIndex(subdomain, view.columns()[col.index()]) != Partition::noIndex)
                    pattern.add(row.index(), partition->localIndex(subdomain, view.columns()[col.index()]));
        pattern.exportIdx(squareBlock_);
        squareBlock_ = 0.0;
    }

    using TimeLoop = typename Assembler::TimeLoop;

    void setTimeLoop(std::shared_ptr<const TimeLoop> timeLoop) { assembler_.setTimeLoop(timeLoop); }
    void setPreviousSolution(const SolutionVector& prevSol) { assembler_.setPreviousSolution(prevSol); }
    void setMaxIterations(std::size_t n) { maxIterations_ = n; }
    void setResidualTolerance(double tol) { residualTolerance_ = tol; }
    void setMaxLineSearchTrials(std::size_t n) { maxLineSearchTrials_ = n; }

    /*!
     * \brief Bound the iterate after an update, mirroring Newton's chopped update
     *
     * Called for every degree of freedom with the newly computed primary variables, the ones they came
     * from, and the degree of freedom's index, so the model can clamp what it knows to be bounded. The
     * index is what makes a spatially varying bound expressible: the pressure range admitted by a
     * saturation-capillary-pressure relation depends on the material at that location, so a bound stated
     * only in terms of the primary variables cannot represent it — a saturation to \f$[0,1]\f$, a
     * pressure to its tabulated range, a per-iteration change to a maximum. Repairing the iterate is
     * preferable to rejecting the step: a single cell that overshoots would otherwise force the whole
     * subdomain's step to be damped, repeatedly, until the local solve gives up.
     */
    void setUpdateBounding(std::function<void(PrimaryVariables&, const PrimaryVariables&, std::size_t)> bound)
    { boundUpdate_ = std::move(bound); }

    /*!
     * \brief Declare which primary variables the model considers admissible
     *
     * Finiteness is always required. A model whose variables carry a range — a saturation, a mole
     * fraction — should say so here, because leaving that range does not necessarily raise an error: it
     * can instead produce a degenerate relative permeability and a singular local block, which surfaces
     * much later and looks like a solver failure rather than an out-of-range state.
     */
    void setPhysicalityCheck(std::function<bool(const PrimaryVariables&)> check)
    { isPhysicalDof_ = std::move(check); }

    //! the switch is owned by the outer solver so that all subdomains share its per-dof bookkeeping
    void setPrimaryVariableSwitch(SubdomainPrimaryVariableSwitch<GridVariables>* pvSwitch)
    { pvSwitch_ = pvSwitch; }

    /*!
     * \brief Solve the subdomain problem in place
     * \param v a global solution vector; only the subdomain's own degrees of freedom are modified
     */
    LocalSolveResult solve(SolutionVector& v)
    {
        LocalSolveResult result;
        const auto& dofs = partition_->dofs(subdomain_);

        LocalVector update(dofs.size());
        LocalVector rhs(dofs.size());
        bool switchedSinceLastCheck = false;

        for (std::size_t k = 0; k <= maxIterations_; ++k)
        {
            try {
                updateSubdomainVariables(*gridGeometry_, *gridVariables_, v, updateRange_);
                assembler_.assembleJacobianAndResidual(v);
            }
            catch (const NumericalProblem&) { return result; }

            // kept in lockstep with the assembly so that a converged solve cannot leave a stale block
            extractSquareBlock_();

            result.residualNorm = residualNorm_();
            if (!std::isfinite(result.residualNorm))
                return result;

            // A degree of freedom that has just changed phase state is described by different primary
            // variables than the residual just assembled was built from one iteration ago, so one more
            // pass is required before the state may be certified as converged.
            if (result.residualNorm < residualTolerance_ && !switchedSinceLastCheck)
            {
                result.converged = true;
                result.iterations = k;
                return result;
            }
            switchedSinceLastCheck = false;

            if (k == maxIterations_)
            {
                result.iterations = k;
                return result;
            }

            for (std::size_t r = 0; r < dofs.size(); ++r)
                for (int eq = 0; eq < int(LocalVector::block_type::dimension); ++eq)
                    rhs[r][eq] = assembler_.residual().vector()[r][eq];

            update = 0.0;
            try {
                Dune::UMFPack<SquareMatrix> solver(squareBlock_, 0);
                Dune::InverseOperatorResult solverResult;
                solver.apply(update, rhs, solverResult);
                if (!solverResult.converged)
                    return result;
            }
            catch (const Dune::Exception&) { return result; }

            // A local problem that is driven far from its previous state by a frozen exterior can
            // overshoot badly on a full step. Backtracking on the local residual keeps such a
            // subdomain converging instead of failing and forcing the whole outer step to be abandoned.
            double lambda = 1.0;
            const auto previousNorm = result.residualNorm;
            SolutionVector trial(v);
            bool accepted = false;
            for (std::size_t trialCount = 0; trialCount <= maxLineSearchTrials_; ++trialCount)
            {
                trial = v;
                for (std::size_t r = 0; r < dofs.size(); ++r)
                {
                    auto& priVars = trial[dofs[r]];
                    for (int eq = 0; eq < int(LocalVector::block_type::dimension); ++eq)
                        priVars[eq] -= lambda*update[r][eq];

                    if (boundUpdate_)
                        boundUpdate_(priVars, v[dofs[r]], dofs[r]);
                }

                if (!isPhysical_(trial))
                {
                    lambda *= 0.5;
                    continue;
                }

                double trialNorm = 0.0;
                try {
                    updateSubdomainVariables(*gridGeometry_, *gridVariables_, trial, updateRange_);
                    assembler_.assembleResidual(trial);
                    trialNorm = residualNorm_();
                }
                catch (const NumericalProblem&)
                {
                    lambda *= 0.5;
                    continue;
                }

                // the last trial is accepted on physicality alone: a heavily damped step that does not
                // yet reduce the residual is still progress, whereas an out-of-range one is not
                if (std::isfinite(trialNorm) && (trialNorm < previousNorm || trialCount == maxLineSearchTrials_))
                {
                    accepted = true;
                    break;
                }

                lambda *= 0.5;
            }

            // no admissible step exists along this direction; reporting failure hands the sweep to the
            // outer solver's fallback rather than continuing from a state outside the model's range
            if (!accepted)
            {
                result.iterations = k;
                return result;
            }

            v = trial;

            // only the owner may switch a degree of freedom, so the range is the core rather than the
            // subdomain: an overlap degree of freedom belongs to a different subdomain's core
            if (pvSwitch_ && pvSwitch_->active())
            {
                updateSubdomainVariables(*gridGeometry_, *gridVariables_, v, updateRange_);
                if (pvSwitch_->applyToSubdomain(v, *gridVariables_, assembler_.problem(),
                                                *gridGeometry_, partition_->core(subdomain_)))
                {
                    switchedSinceLastCheck = true;
                    result.switched += pvSwitch_->numSwitched();
                }
            }
        }

        return result;
    }

    //! valid after a converged solve: the Jacobian assembled at the converged iterate
    const Assembler& assembler() const { return assembler_; }
    const SquareMatrix& squareBlock() const { return squareBlock_; }

private:
    double residualNorm_() const
    {
        double norm = 0.0;
        for (const auto& block : assembler_.residual().vector())
            norm = std::max(norm, static_cast<double>(block.infinity_norm()));
        return norm;
    }

    bool isPhysical_(const SolutionVector& v) const
    {
        for (const auto dof : partition_->dofs(subdomain_))
        {
            const auto& priVars = v[dof];
            for (std::size_t eq = 0; eq < priVars.size(); ++eq)
                if (!std::isfinite(priVars[eq]))
                    return false;

            if (isPhysicalDof_ && !isPhysicalDof_(priVars))
                return false;
        }
        return true;
    }

    void extractSquareBlock_()
    {
        squareBlock_ = 0.0;
        const auto& view = assembler_.jacobian();
        for (auto row = view.matrix().begin(); row != view.matrix().end(); ++row)
            for (auto col = row->begin(); col != row->end(); ++col)
            {
                const auto localCol = partition_->localIndex(subdomain_, view.columns()[col.index()]);
                if (localCol != Partition::noIndex)
                    squareBlock_[row.index()][localCol] = *col;
            }
    }

    Assembler assembler_;
    std::shared_ptr<GridVariables> gridVariables_;
    std::shared_ptr<const GridGeometry> gridGeometry_;
    std::shared_ptr<const Partition> partition_;
    SubdomainIndex subdomain_;
    std::vector<Index> updateRange_;
    SubdomainPrimaryVariableSwitch<GridVariables>* pvSwitch_ = nullptr;
    std::function<bool(const PrimaryVariables&)> isPhysicalDof_;
    std::function<void(PrimaryVariables&, const PrimaryVariables&, std::size_t)> boundUpdate_;

    SquareMatrix squareBlock_;
    std::size_t maxIterations_ = 20;
    std::size_t maxLineSearchTrials_ = 6;
    double residualTolerance_ = 1e-12;
};

} // end namespace Dumux::NonlinearPreconditioning

#endif
