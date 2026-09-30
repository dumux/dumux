// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief A multi-stage (Runge-Kutta) linear system assembler for finite volume schemes
 *        with multiple domains and the new Experimental::GridVariables concept.
 *
 * This class combines multi-stage time integration (Butcher tableau) with the
 * Experimental::MultiDomainAssembler infrastructure (which uses Concept::GridVariables).
 * It is the transient counterpart of Experimental::MultiDomainAssembler.
 *
 * Usage mirrors MultiStageMultiDomainFVAssembler but uses the new GridVariables API:
 * \code
 *   using Assembler = Experimental::MultiStageMultiDomainAssembler<MDTraits, CM, DiffMethod::numeric>;
 *   auto assembler = std::make_shared<Assembler>(problems, gridGeometries, gridVariables,
 *                                                couplingManager, timeSteppingMethod, prevSol);
 *   // in the time loop:
 *   Experimental::MultiStageTimeStepper timeStepper(nonLinearSolver, timeSteppingMethod);
 *   timeStepper.step(x, t, dt);
 * \endcode
 */
#ifndef DUMUX_MULTIDOMAIN_MULTISTAGE_MULTIDOMAIN_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_MULTISTAGE_MULTIDOMAIN_ASSEMBLER_HH

#include <vector>
#include <memory>
#include <type_traits>
#include <tuple>

#include <dune/common/hybridutilities.hh>
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/utility.hh>
#include <dumux/common/gridcapabilities.hh>
#include <dumux/discretization/method.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/jacobianpattern.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/parallel/multithreading.hh>

#include <dumux/multidomain/couplingjacobianpattern.hh>
#include <dumux/multidomain/assemblerview.hh>
#include <dumux/multidomain/subdomaincclocalassembler.hh>
#include <dumux/multidomain/subdomaincvfelocalassembler_.hh>

#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>

// CouplingManagerSupportsMultithreadedAssembly and Detail::hasSubProblemGlobalConstraints
// are defined in multidomain/fvassembler.hh.
// (The old experimental assembler must NOT be included to avoid redeclaration clashes.)
#include <dumux/multidomain/fvassembler.hh>
#include <dune/common/std/type_traits.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief Multi-stage (Runge-Kutta) multi-domain assembler using the new Concept::GridVariables.
 *
 * \tparam MDTraits  MultiDomainTraits (provides SubDomain<i>::GridVariables etc.)
 * \tparam CMType    CouplingManager type
 * \tparam diffMethod  Differentiation method for Jacobian computation
 */
template<class MDTraits, class CMType, DiffMethod diffMethod>
class MultiStageMultiDomainAssembler
{
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

public:
    using Traits = MDTraits;
    using Scalar = typename MDTraits::Scalar;
    using StageParams = Experimental::MultiStageParams<Scalar>;

    template<std::size_t id>
    using LocalResidual = typename MDTraits::template SubDomain<id>::LocalResidual;

    template<std::size_t id>
    using GridVariables = typename MDTraits::template SubDomain<id>::GridVariables;

    template<std::size_t id>
    using GridGeometry = typename MDTraits::template SubDomain<id>::GridGeometry;

    template<std::size_t id>
    using Problem = typename MDTraits::template SubDomain<id>::Problem;

    using JacobianMatrix = typename MDTraits::JacobianMatrix;
    using SolutionVector = typename MDTraits::SolutionVector;
    using ResidualType = typename MDTraits::ResidualVector;
    using CouplingManager = CMType;

private:
    using ProblemTuple = typename MDTraits::template TupleOfSharedPtrConst<Problem>;
    using GridGeometryTuple = typename MDTraits::template TupleOfSharedPtrConst<GridGeometry>;
    using GridVariablesTuple = typename MDTraits::template TupleOfSharedPtr<GridVariables>;
    using ThisType = MultiStageMultiDomainAssembler<MDTraits, CouplingManager, diffMethod>;

    template<std::size_t id>
    using SubDomainAssemblerView = MultiDomainAssemblerSubDomainView<ThisType, id>;

    template<class DiscretizationMethod, std::size_t id>
    struct SubDomainAssemblerType;

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethods::CCTpfa, id>
    {
        using type = SubDomainCCLocalAssembler<id, SubDomainTypeTag<id>, SubDomainAssemblerView<id>, diffMethod, true>;
    };

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethods::CCMpfa, id>
    {
        using type = SubDomainCCLocalAssembler<id, SubDomainTypeTag<id>, SubDomainAssemblerView<id>, diffMethod, true>;
    };

    template<std::size_t id, class DM>
    struct SubDomainAssemblerType<DiscretizationMethods::CVFE<DM>, id>
    {
        using type = Dumux::Experimental::SubDomainCVFELocalAssembler<id, SubDomainTypeTag<id>, SubDomainAssemblerView<id>, diffMethod, true>;
    };

    template<std::size_t id>
    using SubDomainAssembler = typename SubDomainAssemblerType<typename GridGeometry<id>::DiscretizationMethod, id>::type;

public:

    /*!
     * \brief Constructor for instationary multi-stage assembly.
     *
     * \param problem       Tuple of shared_ptr<const Problem> per subdomain
     * \param gridGeometry  Tuple of shared_ptr<const GridGeometry> per subdomain
     * \param gridVariables Tuple of shared_ptr<GridVariables> per subdomain
     * \param couplingManager The coupling manager
     * \param msMethod      The multi-stage time stepping method (Butcher tableau)
     * \param prevSol       The solution at the previous time step
     */
    MultiStageMultiDomainAssembler(ProblemTuple problem,
                                    GridGeometryTuple gridGeometry,
                                    GridVariablesTuple gridVariables,
                                    std::shared_ptr<CouplingManager> couplingManager,
                                    std::shared_ptr<const Experimental::MultiStageMethod<Scalar>> msMethod,
                                    const SolutionVector& prevSol)
    : couplingManager_(couplingManager)
    , timeSteppingMethod_(msMethod)
    , problemTuple_(std::move(problem))
    , gridGeometryTuple_(std::move(gridGeometry))
    , gridVariablesTuple_(std::move(gridVariables))
    , prevSol_(&prevSol)
    {
        std::cout << "Instantiated multi-stage multi-domain assembler." << std::endl;

        enableMultithreading_ = CouplingManagerSupportsMultithreadedAssembly<CouplingManager>::value
            && Grid::Capabilities::allGridsSupportsMultithreading(gridGeometryTuple_)
            && !Multithreading::isSerial()
            && getParam<bool>("Assembly.Multithreading", true);

        maybeComputeColors_();
    }

    /*!
     * \brief Assembles the Jacobian and residual for the current Newton iteration
     *        of the current time integration stage.
     */
    void assembleJacobianAndResidual(const SolutionVector& curSol)
    {
        resetJacobian_();
        resetResidual_();
        spatialOperatorEvaluations_.back() = 0.0;
        temporalOperatorEvaluations_.back() = 0.0;

        if (stageParams_->size() != spatialOperatorEvaluations_.size())
            DUNE_THROW(Dune::InvalidStateException, "Wrong number of stage residuals. Call prepareStage first.");

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
        {
            auto& jacRow = (*jacobian_)[domainId];
            auto& spatial = spatialOperatorEvaluations_.back()[domainId];
            auto& temporal = temporalOperatorEvaluations_.back()[domainId];

            assemble_(domainId, [&](const auto& element)
            {
                MultiDomainAssemblerSubDomainView view{*this, domainId};
                SubDomainAssembler<domainId()> subDomainAssembler(view, element, curSol, *couplingManager_);
                subDomainAssembler.assembleJacobianAndResidual(
                    jacRow, (*residual_)[domainId],
                    gridVariablesTuple_,
                    *stageParams_, temporal, spatial,
                    constrainedDofs_[domainId]
                );
            });

            // add contributions from all previous stages
            auto constantResidualComponent = (*residual_)[domainId];
            constantResidualComponent = 0.0;
            for (std::size_t k = 0; k < stageParams_->size()-1; ++k)
            {
                if (!stageParams_->skipTemporal(k))
                    constantResidualComponent.axpy(stageParams_->temporalWeight(k), temporalOperatorEvaluations_[k][domainId]);
                if (!stageParams_->skipSpatial(k))
                    constantResidualComponent.axpy(stageParams_->spatialWeight(k), spatialOperatorEvaluations_[k][domainId]);
            }

            // apply masking for constrained dofs and add to residual
            for (std::size_t i = 0; i < constantResidualComponent.size(); ++i)
                for (std::size_t ii = 0; ii < constantResidualComponent[i].size(); ++ii)
                    (*residual_)[domainId][i][ii] += constrainedDofs_[domainId][i][ii] > 0.5 ? 0.0 : constantResidualComponent[i][ii];

            // enforce Dirichlet constraints (overwrites Jacobian rows and residual for constrained DOFs)
            enforceProblemConstraints_(domainId, (*jacobian_)[domainId], (*residual_)[domainId],
                                       gridGeometry(domainId), curSol[domainId]);
        });
    }

    //! assembleResidual is not supported (only Jacobian+residual combined)
    void assembleResidual(const SolutionVector&)
    { DUNE_THROW(Dune::NotImplemented, "assembleResidual without Jacobian not implemented for multi-stage assembly"); }

    //! Set up internal Jacobian and residual storage
    void setLinearSystem()
    {
        jacobian_ = std::make_shared<JacobianMatrix>();
        residual_ = std::make_shared<ResidualType>();

        setJacobianBuildMode_(*jacobian_);
        setJacobianPattern_(*jacobian_);
        setResidualSize_(*residual_);
    }

    //! Update all grid variables and the coupling manager from a new solution.
    void updateGridVariables(const SolutionVector& curSol)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(gridVariablesTuple_)), [&](const auto domainId)
        { this->gridVariables(domainId).update(curSol[domainId]); });

        // Keep coupling manager in sync so cross-domain queries (e.g. J from displacement)
        // always reflect the current Newton iterate.
        couplingManager_->updateSolution(curSol);
    }

    //! Reset grid variables to a previous state (e.g., on Newton failure)
    void resetTimeStep(const SolutionVector& curSol)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(gridVariablesTuple_)), [&](const auto domainId)
        { this->gridVariables(domainId).resetTimeStep(curSol[domainId]); });

        this->clearStages();
    }

    //! Prepare the assembler for stage \p curStage+1 of the time integration.
    template<class StageParamsPtr>
    void prepareStage(SolutionVector& x, StageParamsPtr params)
    {
        stageParams_ = std::move(params);
        const auto curStage = stageParams_->size() - 1;

        // At stage 1: also assemble the stage-0 (previous time level) residuals.
        if (curStage == 1)
        {
            using namespace Dune::Hybrid;
            forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
            {
                setProblemTime_(*std::get<domainId>(problemTuple_), stageParams_->timeAtStage(0));
            });

            resetResidual_();
            spatialOperatorEvaluations_.push_back(*residual_);
            temporalOperatorEvaluations_.push_back(*residual_);

            // assemble stage 0 using x = prevSol
            forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
            {
                auto& spatial = spatialOperatorEvaluations_.back()[domainId];
                auto& temporal = temporalOperatorEvaluations_.back()[domainId];
                assemble_(domainId, [&](const auto& element)
                {
                    MultiDomainAssemblerSubDomainView view{*this, domainId};
                    SubDomainAssembler<domainId()> subDomainAssembler(view, element, x, *couplingManager_);
                    subDomainAssembler.assembleCurrentResidual(temporal, spatial);
                });
            });
        }

        // update time on all problems for the current stage
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
        {
            setProblemTime_(*std::get<domainId>(problemTuple_), stageParams_->timeAtStage(curStage));
        });

        resetResidual_();
        spatialOperatorEvaluations_.push_back(*residual_);
        temporalOperatorEvaluations_.push_back(*residual_);
    }

    //! Clear all stored stage evaluations (call at the start of each time step)
    void clearStages()
    {
        spatialOperatorEvaluations_.clear();
        temporalOperatorEvaluations_.clear();
        stageParams_.reset();
    }

    //! Whether this assembler is for a stationary problem (always false here)
    bool isStationaryProblem() const { return false; }
    //! Whether the time integration is implicit (always true for multi-stage)
    bool isImplicit() const { return timeSteppingMethod_->implicit(); }

    //! Number of dofs of subdomain i
    template<std::size_t i>
    std::size_t numDofs(Dune::index_constant<i> domainId) const
    { return std::get<domainId>(gridGeometryTuple_)->numDofs(); }

    template<std::size_t i> const auto& problem(Dune::index_constant<i> domainId) const
    { return *std::get<domainId>(problemTuple_); }

    template<std::size_t i> const auto& gridGeometry(Dune::index_constant<i> domainId) const
    { return *std::get<domainId>(gridGeometryTuple_); }

    template<std::size_t i> const auto& gridView(Dune::index_constant<i> domainId) const
    { return gridGeometry(domainId).gridView(); }

    template<std::size_t i> GridVariables<i>& gridVariables(Dune::index_constant<i> domainId)
    { return *std::get<domainId>(gridVariablesTuple_); }

    template<std::size_t i> const GridVariables<i>& gridVariables(Dune::index_constant<i> domainId) const
    { return *std::get<domainId>(gridVariablesTuple_); }

    const CouplingManager& couplingManager() const { return *couplingManager_; }

    JacobianMatrix& jacobian() { return *jacobian_; }
    ResidualType& residual() { return *residual_; }

    const SolutionVector& prevSol() const { return *prevSol_; }

    //! Update the pointer to the previous-time-step solution (call after each accepted step).
    void setPreviousSolution(const SolutionVector& sol)
    { prevSol_ = &sol; }

    //! Returns the current stage temporal and spatial weights (for use in evalCouplingResidual).
    std::pair<Scalar, Scalar> currentStageWeights() const
    {
        if (!stageParams_)
            DUNE_THROW(Dune::InvalidStateException, "No stage params set. Call prepareStage first.");
        const auto k = stageParams_->size() - 1;
        return {stageParams_->temporalWeight(k), stageParams_->spatialWeight(k)};
    }

    template<std::size_t i>
    LocalResidual<i> localResidual(Dune::index_constant<i> domainId) const
    { return LocalResidual<i>(std::get<domainId>(problemTuple_).get(), nullptr); }

    template<std::size_t i>
    auto& curSol(Dune::index_constant<i> id) { return (*curSol_)[id]; }
    template<std::size_t i>
    const auto& curSol(Dune::index_constant<i> id) const { return (*curSol_)[id]; }

    void computeColorsForAssembly()
    {
        if constexpr (CouplingManagerSupportsMultithreadedAssembly<CouplingManager>::value)
            couplingManager_->computeColorsForAssembly();
    }

    template<std::size_t i, class AssembleElementFunc>
    void assembleMultithreaded(Dune::index_constant<i>, AssembleElementFunc&& f) const
    {
        if constexpr (CouplingManagerSupportsMultithreadedAssembly<CouplingManager>::value)
            couplingManager_->assembleMultithreaded(Dune::index_constant<i>{}, std::forward<AssembleElementFunc>(f));
    }

protected:
    std::shared_ptr<CouplingManager> couplingManager_;

private:
    void setJacobianBuildMode_(JacobianMatrix& jac) const
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto i)
        {
            forEach(jac[i], [&](auto& jacBlock)
            {
                using BlockType = std::decay_t<decltype(jacBlock)>;
                if (jacBlock.buildMode() == BlockType::BuildMode::unknown)
                    jacBlock.setBuildMode(BlockType::BuildMode::random);
                else if (jacBlock.buildMode() != BlockType::BuildMode::random)
                    DUNE_THROW(Dune::NotImplemented, "Only BCRS matrices with random build mode are supported");
            });
        });
    }

    void setJacobianPattern_(JacobianMatrix& jac) const
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainI)
        {
            forEach(integralRange(Dune::Hybrid::size(jac[domainI])), [&](const auto domainJ)
            {
                const auto pattern = getJacobianPattern_(domainI, domainJ);
                pattern.exportIdx(jac[domainI][domainJ]);
            });
        });
    }

    void setResidualSize_(ResidualType& res) const
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(res)), [&](const auto domainId)
        { res[domainId].resize(this->numDofs(domainId)); });
    }

    void resetResidual_()
    {
        if (!residual_)
        {
            residual_ = std::make_shared<ResidualType>();
            setResidualSize_(*residual_);
        }

        setResidualSize_(constrainedDofs_);
        (*residual_) = 0.0;
        constrainedDofs_ = 0.0;
    }

    void resetJacobian_()
    {
        if (!jacobian_)
        {
            jacobian_ = std::make_shared<JacobianMatrix>();
            setJacobianBuildMode_(*jacobian_);
            setJacobianPattern_(*jacobian_);
        }
        (*jacobian_) = 0.0;
    }

    void maybeComputeColors_()
    {
        if constexpr (CouplingManagerSupportsMultithreadedAssembly<CouplingManager>::value)
            if (enableMultithreading_)
                couplingManager_->computeColorsForAssembly();
    }

    template<std::size_t i, class AssembleElementFunc>
    void assemble_(Dune::index_constant<i> domainId, AssembleElementFunc&& assembleElement) const
    {
        bool succeeded = false;
        try
        {
            if constexpr (CouplingManagerSupportsMultithreadedAssembly<CouplingManager>::value)
            {
                if (enableMultithreading_)
                {
                    couplingManager_->assembleMultithreaded(
                        domainId, std::forward<AssembleElementFunc>(assembleElement)
                    );
                    return;
                }
            }

            for (const auto& element : elements(gridView(domainId)))
                assembleElement(element);

            succeeded = true;
        }
        catch (NumericalProblem& e)
        {
            std::cout << "rank " << gridView(domainId).comm().rank()
                      << " caught an exception while assembling: " << e.what() << "\n";
            succeeded = false;
        }

        if (gridView(domainId).comm().size() > 1)
            succeeded = gridView(domainId).comm().min(succeeded);

        if (!succeeded)
            DUNE_THROW(NumericalProblem, "A process did not succeed in linearizing the system");
    }

    //! Enforce Dirichlet constraints from problem.constraints() into Jacobian and residual.
    //! \param jacRow  The full Jacobian row for domainI (MultiTypeBlockVector of Jacobian blocks)
    //! \param res     The per-domain residual vector
    template<std::size_t i, class JacRow, class Res, class GG, class Sol>
    void enforceProblemConstraints_(Dune::index_constant<i> domainI,
                                     JacRow& jacRow, Res& res,
                                     const GG& /*gridGeometry*/, const Sol& curSol) const
    {
        if constexpr (Dumux::Detail::hasSubProblemGlobalConstraints<Problem<domainI>>())
        {
            auto& jac = jacRow[domainI]; // diagonal BCRSMatrix block for this domain

            auto applyDirichletConstraint = [&](const auto& dofIdx, const auto& values,
                                                 const auto eqIdx, const auto pvIdx)
            {
                res[dofIdx][eqIdx] = curSol[dofIdx][pvIdx] - values[pvIdx];

                // zero the diagonal row for this equation
                auto& row = jac[dofIdx];
                for (auto col = row.begin(); col != row.end(); ++col)
                    row[col.index()][eqIdx] = 0.0;
                jac[dofIdx][dofIdx][eqIdx][pvIdx] = 1.0;

                // also zero the coupling block rows so the full Jacobian row is consistent
                using namespace Dune::Hybrid;
                forEach(makeIncompleteIntegerSequence<JacobianMatrix::N(), domainI>(), [&](const auto couplingId)
                {
                    auto& jacCoupling = jacRow[couplingId];
                    auto& rowCoupling = jacCoupling[dofIdx];
                    for (auto c = rowCoupling.begin(); c != rowCoupling.end(); ++c)
                        rowCoupling[c.index()][eqIdx] = 0.0;
                });
            };

            for (const auto& cd : this->problem(domainI).constraints())
            {
                const auto& info = cd.constraintInfo();
                const auto& vals = cd.values();
                const auto dofIdx = cd.dofIndex();
                for (int eqIdx = 0; eqIdx < info.size(); ++eqIdx)
                    if (info.isConstraintEquation(eqIdx))
                        applyDirichletConstraint(dofIdx, vals, eqIdx, info.eqToPriVarIndex(eqIdx));
            }
        }
    }

    template<std::size_t i, std::size_t j, typename std::enable_if_t<(i==j), int> = 0>
    Dune::MatrixIndexSet getJacobianPattern_(Dune::index_constant<i> domainI,
                                              Dune::index_constant<j>) const
    {
        const auto& gg = gridGeometry(domainI);
        if (timeSteppingMethod_->implicit())
        {
            auto pattern = getJacobianPattern<true>(gg);
            couplingManager_->extendJacobianPattern(domainI, pattern);
            return pattern;
        }
        else
        {
            auto pattern = getJacobianPattern<false>(gg);
            couplingManager_->extendJacobianPattern(domainI, pattern);
            return pattern;
        }
    }

    template<std::size_t i, std::size_t j, typename std::enable_if_t<(i!=j), int> = 0>
    Dune::MatrixIndexSet getJacobianPattern_(Dune::index_constant<i> domainI,
                                              Dune::index_constant<j> domainJ) const
    {
        if (timeSteppingMethod_->implicit())
            return getCouplingJacobianPattern<true>(*couplingManager_,
                domainI, gridGeometry(domainI),
                domainJ, gridGeometry(domainJ));
        else
            return getCouplingJacobianPattern<false>(*couplingManager_,
                domainI, gridGeometry(domainI),
                domainJ, gridGeometry(domainJ));
    }

    template<class P>
    void setProblemTime_(const P& p, const Scalar t)
    { setProblemTimeImpl_(p, t, 0); }

    template<class P>
    auto setProblemTimeImpl_(const P& p, const Scalar t, int) -> decltype(p.setTime(0))
    { p.setTime(t); }

    template<class P>
    void setProblemTimeImpl_(const P& p, const Scalar t, long) {}

    std::shared_ptr<const Experimental::MultiStageMethod<Scalar>> timeSteppingMethod_;
    std::vector<ResidualType> spatialOperatorEvaluations_;
    std::vector<ResidualType> temporalOperatorEvaluations_;
    ResidualType constrainedDofs_;
    std::shared_ptr<const StageParams> stageParams_;

    ProblemTuple problemTuple_;
    GridGeometryTuple gridGeometryTuple_;
    GridVariablesTuple gridVariablesTuple_;

    const SolutionVector* prevSol_ = nullptr;
    const SolutionVector* curSol_ = nullptr;

    std::shared_ptr<JacobianMatrix> jacobian_;
    std::shared_ptr<ResidualType> residual_;

    bool enableMultithreading_ = false;
};

} // end namespace Dumux::Experimental

#endif
