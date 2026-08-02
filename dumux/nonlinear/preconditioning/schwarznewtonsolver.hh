// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Nonlinear
 * \brief Additive Schwarz preconditioned exact Newton for non-overlapping subdomains
 */
#ifndef DUMUX_NONLINEAR_PRECONDITIONING_SCHWARZ_NEWTON_SOLVER_HH
#define DUMUX_NONLINEAR_PRECONDITIONING_SCHWARZ_NEWTON_SOLVER_HH

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/umfpack.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/linear/dunevectors.hh>

#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/nonlinear/preconditioning/localsolver.hh>
#include <dumux/nonlinear/preconditioning/partition.hh>
#include <dumux/nonlinear/preconditioning/partitioner.hh>
#include <dumux/nonlinear/preconditioning/pvswitch.hh>

namespace Dumux::NonlinearPreconditioning {

struct SchwarzNewtonReport
{
    bool converged = false;
    std::size_t outerIterations = 0;
    std::size_t globalStages = 0;      //!< outer iterations that needed a global solve
    std::size_t localIterations = 0;
    std::size_t failedLocalSolves = 0;
    std::size_t newtonFallbackSteps = 0;  //!< outer steps taken on the true system after a local failure
    std::size_t switchedDofs = 0;         //!< degrees of freedom that changed phase state
    std::size_t rejectedOuterSteps = 0;   //!< outer steps abandoned because no admissible damping existed
    double residualNorm = 0.0;
};

/*!
 * \ingroup Nonlinear
 * \brief Solves \f$ F(u) = 0 \f$ by nonlinear preconditioning with non-overlapping subdomain solves
 *
 * Each outer iteration has a local and a global stage. The local stage solves every subdomain
 * independently with the other subdomains' unknowns frozen, producing the corrected iterates
 * \f$ u^{(i)} \f$ and, as a by-product, the Jacobian rows \f$ R_i J(u^{(i)}) \f$.
 *
 * The global stage would solve \f$ \partial\tilde{F}/\partial u\,\Delta u = \tilde{F}(u) \f$, but that
 * operator is dense. For a non-overlapping decomposition the block diagonal
 * \f$ D = \mathrm{blockdiag}(R_i J(u^{(i)}) P_i) \f$ is exact, and multiplying through by it gives the
 * equivalent system
 * \f[ -\frac{\partial R}{\partial u}\,\Delta u = D\,\tilde{F}(u), \f]
 * whose matrix is an ordinary sparse global Jacobian, assembled row-block by row-block from the local
 * stage. This is solvable with the same linear solvers and preconditioners as an unpreconditioned Newton
 * step, which is why it is preferred here over the matrix-free form.
 *
 * The global residual is checked after the local stage. Besides being the guard against accepting a step
 * whose frozen regions are not actually converged, this frequently lets the global stage be skipped
 * altogether.
 */
template<class TypeTag, DiffMethod diffMethod = DiffMethod::numeric>
class SchwarzNewtonSolver
{
    using LocalSolver = SubdomainSolver<TypeTag, diffMethod>;
    using GlobalAssembler = FVAssembler<TypeTag, diffMethod>;

public:
    using SolutionVector = typename LocalSolver::SolutionVector;
    using PrimaryVariables = typename LocalSolver::PrimaryVariables;
    using Problem = typename LocalSolver::Problem;
    using GridGeometry = typename LocalSolver::GridGeometry;
    using GridVariables = typename LocalSolver::GridVariables;
    using Scalar = typename GlobalAssembler::Scalar;
    using TimeLoop = TimeLoopBase<Scalar>;
    using GlobalMatrix = typename GlobalAssembler::JacobianMatrix;
    using ResidualVector = typename GlobalAssembler::ResidualType;
    using Index = Partition::Index;

    SchwarzNewtonSolver(std::shared_ptr<const Problem> problem,
                std::shared_ptr<const GridGeometry> gridGeometry,
                std::shared_ptr<GridVariables> gridVariables,
                std::shared_ptr<const Partition> partition,
                const std::vector<std::vector<Index>>& influence,
                const std::string& paramGroup = "")
    : partition_(partition)
    , gridVariables_(gridVariables)
    , problem_(problem)
    , gridGeometry_(gridGeometry)
    , pvSwitch_(paramGroup)
    {
        if (partition->overlapWidth() != 0)
            DUNE_THROW(Dune::InvalidStateException,
                       "The ASPEN reformulation requires a non-overlapping partition, but the overlap "
                       "width is " << partition->overlapWidth());

        // The switch is deliberately not wired into the local solves by default; see
        // setSwitchInLocalSolves for why the invocation cadence matters.
        for (SubdomainIndex i = 0; i < partition->numSubdomains(); ++i)
            solvers_.emplace_back(std::make_unique<LocalSolver>(
                problem, gridGeometry, gridVariables, partition, influence, i));

        globalAssembler_ = std::make_shared<GlobalAssembler>(problem, gridGeometry, gridVariables);

        Dune::MatrixIndexSet pattern(partition->numDofs(), partition->numDofs());
        for (Index col = 0; col < static_cast<Index>(influence.size()); ++col)
        {
            pattern.add(col, col);
            for (const auto row : influence[col])
                pattern.add(row, col);
        }
        pattern.exportIdx(matrix_);
    }

    void setTimeLoop(std::shared_ptr<const TimeLoop> timeLoop)
    {
        timeLoop_ = timeLoop;
        for (auto& s : solvers_)
            s->setTimeLoop(timeLoop);
    }

    void setPreviousSolution(const SolutionVector& prevSol)
    {
        prevSol_ = &prevSol;
        for (auto& s : solvers_)
            s->setPreviousSolution(prevSol);
    }

    void setMaxIterations(std::size_t n) { maxIterations_ = n; }
    void setMaxLineSearchTrials(std::size_t n) { maxLineSearchTrials_ = n; }
    /*!
     * \brief Absolute floor on the true global residual, in the residual's own units
     * \note Disabled by default. An absolute threshold is only meaningful once the residual's scale is
     *       known, and that scale varies by many orders of magnitude between problems: a fine grid with
     *       a small source term produces residuals that sit below any fixed default from the outset,
     *       which would let the gate accept a step on which nothing had happened.
     */
    void setResidualTolerance(double tol) { residualTolerance_ = tol; }

    //! Required reduction of the true global residual relative to its value at the start of the solve
    void setResidualReduction(double reduction) { residualReduction_ = reduction; }

    //! Relative movement of the solution below which a sweep counts as settled, as Newton uses
    void setMaxRelativeShift(double shift) { maxRelativeShift_ = shift; }
    void setToleratedPriVarSwitches(std::size_t n) { toleratedSwitches_ = n; }

    //! bound the iterate after an update, in both stages; see SubdomainSolver::setUpdateBounding
    void setUpdateBounding(std::function<void(PrimaryVariables&, const PrimaryVariables&, std::size_t)> bound)
    {
        boundUpdate_ = bound;
        for (auto& s : solvers_)
            s->setUpdateBounding(bound);
    }

    //! declare which primary variables the model considers admissible; see SubdomainSolver
    void setPhysicalityCheck(std::function<bool(const PrimaryVariables&)> check)
    {
        isPhysicalDof_ = check;
        for (auto& s : solvers_)
            s->setPhysicalityCheck(check);
    }

    /*!
     * \brief Also apply the phase switch inside every local Newton iteration
     *
     * Off by default, and enabling it is not free. A switch implementation records whether a degree of
     * freedom crossed its bare threshold on *every* call, and uses that record to widen the threshold on
     * the next one — hysteresis sized for Newton's cadence of roughly one call per iteration. A local
     * stage calls it once per local iteration per subdomain, which is orders of magnitude more often, so
     * a cell merely hovering near the threshold latches the widened bar and stops being able to cross
     * it. The visible symptom is a phase front that stops advancing while the solver reports clean
     * convergence throughout.
     *
     * Applying the switch once per outer sweep to the combined solution, which the solver always does,
     * matches Newton's cadence and does not have this interaction.
     */
    void setSwitchInLocalSolves(bool val)
    {
        for (auto& s : solvers_)
            s->setPrimaryVariableSwitch(val ? &pvSwitch_ : nullptr);
    }
    void setVerbosity(int v) { verbosity_ = v; }
    void setLocalMaxIterations(std::size_t n)
    {
        for (auto& s : solvers_)
            s->setMaxIterations(n);
    }
    void setLocalResidualTolerance(double tol)
    {
        for (auto& s : solvers_)
            s->setResidualTolerance(tol);
    }

    SchwarzNewtonReport solve(SolutionVector& u)
    {
        SchwarzNewtonReport report;

        globalAssembler_->setLinearSystem();
        if (timeLoop_)
        {
            globalAssembler_->setTimeLoop(timeLoop_);
            globalAssembler_->setPreviousSolution(*prevSol_);
        }

        pvSwitch_.initialize(u, *gridVariables_, *problem_, *gridGeometry_);

        // Measured at the incoming iterate, before the first local stage has changed anything. Taking
        // it after that stage would make the target a reduction of an already-reduced residual, which
        // for a single subdomain — where the local solve is the whole problem — is unreachable.
        try {
            gridVariables_->update(u);
            globalAssembler_->assembleResidual(u);
            initialResidualNorm_ = globalResidualNorm_();
        }
        catch (const NumericalProblem&) { initialResidualNorm_ = 0.0; }

        SolutionVector preconditioned(u), combined(u);
        ResidualVector update(u.size());
        scratch_ = u;

        for (std::size_t k = 0; k < maxIterations_; ++k)
        {
            report.outerIterations = k + 1;

            // local stage: every subdomain is solved with the others frozen
            preconditioned = u;
            preconditioned = 0.0;
            combined = u;
            scratch_ = u;
            const auto failuresBefore = report.failedLocalSolves;
            std::size_t switchedThisSweep = 0;
            for (SubdomainIndex i = 0; i < partition_->numSubdomains(); ++i)
            {
                // the previous subdomain only touched its own degrees of freedom, so the scratch vector
                // is restored in O(n_i) rather than rebuilt from u at O(n) per subdomain
                if (i > 0)
                    for (const auto dof : partition_->dofs(i - 1))
                        scratch_[dof] = u[dof];

                auto& v = scratch_;
                const auto local = solvers_[i]->solve(v);
                report.localIterations += local.iterations;
                switchedThisSweep += local.switched;

                // A subdomain that fails to converge contributes no correction this sweep. That is a
                // legitimate zero contribution rather than an error: the outer iteration still has the
                // global residual to converge, and the true residual gate below decides whether the
                // step may be accepted. Failing the whole solve here would throw away the progress
                // every other subdomain just made.
                if (!local.converged)
                {
                    ++report.failedLocalSolves;
                    continue;
                }

                for (const auto dof : partition_->core(i))
                {
                    combined[dof] = v[dof];
                    preconditioned[dof] = u[dof] - v[dof];
                }
            }

            // The check is made at the combined subdomain solution rather than at the current iterate:
            // that vector is one step of the nonlinear Schwarz fixed point, so accepting it when it is
            // already converged costs nothing and removes a global solve. Every degree of freedom is
            // checked, so a subdomain the local stage left alone cannot hide an unconverged residual
            // behind a converged preconditioned one.
            // Mirror what Newton does at the end of every step: the switch is applied to the candidate
            // iterate before it is judged. The local stage only ever switches degrees of freedom its own
            // subdomain owns and only while its solve is running, so a transition that becomes due once
            // the subdomain solutions are combined is caught nowhere else. Doing it here rather than in
            // one downstream branch also means the accepting path cannot return a solution on which a
            // transition is still pending, which no convergence diagnostic would reveal.
            gridVariables_->update(combined);
            if (pvSwitch_.applyGlobal(combined, *gridVariables_, *problem_, *gridGeometry_))
                switchedThisSweep += pvSwitch_.numSwitched();

            try {
                gridVariables_->update(combined);
                globalAssembler_->assembleResidual(combined);
                report.residualNorm = globalResidualNorm_();
            }
            catch (const NumericalProblem&)
            {
                report.residualNorm = std::numeric_limits<double>::infinity();
            }
            if (verbosity_ > 0)
                std::cout << "    outer " << k << ": global residual = " << report.residualNorm
                          << ", failed local solves = " << report.failedLocalSolves
                          << ", switched this sweep = " << switchedThisSweep << std::endl;
            report.switchedDofs += switchedThisSweep;

            // a sweep that is still resolving phase transitions has not settled, however small the
            // residual happens to look at this particular state
            // Scale-free by construction: the bar is a reduction of this solve's own initial residual,
            // with the absolute floor only applying if one was set. A fixed absolute threshold cannot
            // serve here, because whether it is loose or strict depends entirely on the problem's cell
            // size, source magnitude and time step, and when it is loose this gate accepts on the first
            // sweep with the solution untouched.
            const double bar = std::max(residualTolerance_, residualReduction_*initialResidualNorm_);
            const bool residualConverged = report.residualNorm <= bar
                                           || initialResidualNorm_ <= 0.0;

            // How far a residual can be driven is bounded by floating point, and where that floor lies
            // is a property of the problem: on a fine grid with a small source term only four or five
            // orders may be available, so a reduction target below the floor is unreachable and the
            // solve burns its whole budget failing to meet it. A settled solution is the other, equally
            // valid, statement that there is nothing left to do — it is Newton's own default criterion,
            // and it is immune to that floor. It is accepted only once the residual has genuinely come
            // down, so that a sweep which never moved anything cannot qualify.
            const double shift = Detail::Newton::maxRelativeShift<double>(u, combined);
            const bool settled = shift <= maxRelativeShift_
                                 && report.residualNorm <= stagnationFactor_*initialResidualNorm_;

            // A sweep in which every subdomain failed produced no correction anywhere, so the combined
            // solution is the iterate we started from. Accepting it would report success for a step that
            // did nothing, which is what a residual already below the bar would otherwise allow.
            const bool madeNoProgress = report.failedLocalSolves - failuresBefore
                                        == partition_->numSubdomains();

            if ((residualConverged || settled) && switchedThisSweep <= toleratedSwitches_ && !madeNoProgress)
            {
                u = combined;
                report.converged = true;
                return report;
            }

            // A sweep in which phases appeared or disappeared has done two things: it moved the solution,
            // and it changed which equations some degrees of freedom are governed by. The second is only
            // recorded in the subdomains' own working vectors, so it must be committed by taking the
            // combined solution as the next iterate — otherwise every sweep rediscovers the same
            // transitions and the iteration never terminates. Stepping the preconditioned system instead
            // would also mean trusting a Jacobian across a surface on which the local corrections are not
            // differentiable, which it does not model.
            const bool localFailure = report.failedLocalSolves > failuresBefore;

            // Phase-state changes are committed whenever they occur, not only when they exceed what
            // convergence tolerates: they exist only in the combined solution, and discarding them
            // means rediscovering the same transitions on every subsequent sweep. The tolerance governs
            // when a sweep may be called converged, and nothing else.
            if (switchedThisSweep > 0)
            {
                u = combined;
                gridVariables_->update(u);
                pvSwitch_.applyGlobal(u, *gridVariables_, *problem_, *gridGeometry_);
            }

            // A failed local solve is checked independently of switching rather than being skipped by
            // an early exit from the branch above. When both occur — most sweeps, while a front is
            // advancing through a compositional problem — the fallback must still run, otherwise the
            // stall it exists to prevent is exactly what happens.
            if (switchedThisSweep > 0 && !localFailure)
                continue;

            ++report.globalStages;

            // A subdomain that failed to converge leaves behind a Jacobian row block that was never
            // evaluated at a corrected iterate. Feeding that into the preconditioned system corrupts the
            // step, and because the correction for that subdomain is also zero the iteration can stall
            // exactly rather than fail. When any local solve failed, take an ordinary Newton step on the
            // true system instead: it is the fallback that makes the method no worse than plain Newton.
            if (localFailure)
            {
                ++report.newtonFallbackSteps;
                globalAssembler_->assembleJacobianAndResidual(u);
                matrix_ = globalAssembler_->jacobian();
                rhs_ = globalAssembler_->residual();
            }
            else
                assembleGlobalStage_(preconditioned);

            update = 0.0;
            try {
                Dune::UMFPack<GlobalMatrix> solver(matrix_, 0);
                Dune::InverseOperatorResult solverResult;
                solver.apply(update, rhs_, solverResult);
                if (!solverResult.converged)
                    return report;
            }
            catch (const Dune::Exception&) { return report; }

            // The outer step is a Newton step on the preconditioned residual and is no more globally
            // convergent than any other; without damping a single overshoot far from the solution can
            // leave the iteration wandering for the rest of its budget. Backtracking is measured on the
            // true global residual, which is assembled here anyway for the convergence gate.
            {
                const auto referenceNorm = report.residualNorm;
                double lambda = 1.0;
                SolutionVector trial(u);
                bool accepted = false;
                for (std::size_t trialCount = 0; trialCount <= maxLineSearchTrials_; ++trialCount)
                {
                    trial = u;
                    for (std::size_t i = 0; i < trial.size(); ++i)
                    {
                        auto& priVars = trial[i];
                        for (int eq = 0; eq < int(ResidualVector::block_type::dimension); ++eq)
                            priVars[eq] -= lambda*update[i][eq];

                        if (boundUpdate_)
                            boundUpdate_(priVars, u[i], i);
                    }

                    if (!isPhysical_(trial))
                    {
                        lambda *= 0.5;
                        continue;
                    }

                    double trialNorm = 0.0;
                    try {
                        gridVariables_->update(trial);
                        globalAssembler_->assembleResidual(trial);
                        trialNorm = globalResidualNorm_();
                    }
                    catch (const NumericalProblem&)
                    {
                        lambda *= 0.5;
                        continue;
                    }

                    if (std::isfinite(trialNorm)
                        && (trialNorm < referenceNorm || trialCount == maxLineSearchTrials_))
                    {
                        accepted = true;
                        break;
                    }

                    lambda *= 0.5;
                }

                // the outer iterate is never moved outside the model's range; leaving u where it was
                // costs one wasted sweep, whereas accepting an inadmissible state corrupts every
                // subsequent local solve
                if (accepted)
                    u = trial;
                else
                    report.rejectedOuterSteps++;
            }
            gridVariables_->update(u);
            pvSwitch_.applyGlobal(u, *gridVariables_, *problem_, *gridGeometry_);
        }

        return report;
    }

private:
    bool isPhysical_(const SolutionVector& u) const
    {
        for (std::size_t i = 0; i < u.size(); ++i)
        {
            const auto& priVars = u[i];
            for (std::size_t eq = 0; eq < priVars.size(); ++eq)
                if (!std::isfinite(priVars[eq]))
                    return false;

            if (isPhysicalDof_ && !isPhysicalDof_(priVars))
                return false;
        }
        return true;
    }

    double globalResidualNorm_() const
    {
        double norm = 0.0;
        for (const auto& block : globalAssembler_->residual())
            norm = std::max(norm, static_cast<double>(block.infinity_norm()));
        return norm;
    }

    /*!
     * \brief Build \f$ \partial R/\partial u \f$ and the right-hand side \f$ D\tilde{F}(u) \f$
     *
     * Row block i comes from the local stage, already evaluated at the corrected iterate, so no global
     * reassembly is needed. The right-hand side applies the same subdomain's diagonal block, which is
     * what turns the dense preconditioned operator into this sparse one.
     */
    void assembleGlobalStage_(const SolutionVector& preconditioned)
    {
        matrix_ = 0.0;
        rhs_.resize(preconditioned.size());
        rhs_ = 0.0;

        for (SubdomainIndex i = 0; i < partition_->numSubdomains(); ++i)
        {
            const auto& view = solvers_[i]->assembler().jacobian();
            const auto& columns = view.columns();
            const auto& rows = view.rows();

            for (auto row = view.matrix().begin(); row != view.matrix().end(); ++row)
                for (auto col = row->begin(); col != row->end(); ++col)
                    matrix_[rows[row.index()]][columns[col.index()]] = *col;

            const auto& diagonal = solvers_[i]->squareBlock();
            for (auto row = diagonal.begin(); row != diagonal.end(); ++row)
                for (auto col = row->begin(); col != row->end(); ++col)
                {
                    const auto& src = preconditioned[rows[col.index()]];
                    typename ResidualVector::block_type plain;
                    for (int eq = 0; eq < int(plain.dimension); ++eq)
                        plain[eq] = src[eq];
                    col->umv(plain, rhs_[rows[row.index()]]);
                }
        }
    }

    std::shared_ptr<const Partition> partition_;
    std::shared_ptr<GridVariables> gridVariables_;
    std::shared_ptr<const Problem> problem_;
    std::shared_ptr<const GridGeometry> gridGeometry_;
    SubdomainPrimaryVariableSwitch<GridVariables> pvSwitch_;
    std::vector<std::unique_ptr<LocalSolver>> solvers_;
    std::shared_ptr<GlobalAssembler> globalAssembler_;

    GlobalMatrix matrix_;
    ResidualVector rhs_;
    SolutionVector scratch_;

    std::shared_ptr<const TimeLoop> timeLoop_;
    const SolutionVector* prevSol_ = nullptr;
    std::size_t maxIterations_ = 30;
    std::size_t maxLineSearchTrials_ = 8;
    double residualTolerance_ = 0.0;
    double residualReduction_ = 1e-5;   //!< as Dumux::NewtonSolver's own default
    double maxRelativeShift_ = 1e-8;    //!< as Dumux::NewtonSolver's own default
    double stagnationFactor_ = 1e-2;    //!< residual must have fallen this far before a settled solution counts
    double initialResidualNorm_ = 0.0;
    int verbosity_ = 0;
    std::size_t toleratedSwitches_ = 0;
    std::function<bool(const PrimaryVariables&)> isPhysicalDof_;
    std::function<void(PrimaryVariables&, const PrimaryVariables&, std::size_t)> boundUpdate_;
};

} // end namespace Dumux::NonlinearPreconditioning

#endif
