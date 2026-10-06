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
 * \brief A multi-stage (Runge-Kutta) linear system assembler for multiple domains
 *        discretized with the control-volume finite element interface.
 *
 * This is the multi-stage counterpart of Experimental::MultiDomainAssembler.
 * \code
 *   using Assembler = Experimental::MultiStageMultiDomainAssembler<MDTraits, CM, DiffMethod::numeric>;
 *   auto assembler = std::make_shared<Assembler>(problems, gridDiscretizations, gridVariables,
 *                                                couplingManager, timeSteppingMethod, prevSol);
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
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/utility.hh>
#include <dumux/common/gridcapabilities.hh>
#include <dumux/discretization/method.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/jacobianpattern.hh>
#include <dumux/parallel/multithreading.hh>

#include <dumux/multidomain/couplingjacobianpattern.hh>
#include <dumux/multidomain/assemblerview.hh>
#include <dumux/multidomain/subdomaincvfelocalassembler_.hh>
#include <dumux/multidomain/assemblytraits.hh>

#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief Multi-stage (Runge-Kutta) multi-domain assembler for the control-volume finite element interface
 *
 * \tparam MDTraits MultiDomainTraits
 * \tparam CMType The coupling manager type
 * \tparam diffMethod The differentiation method for the Jacobian
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
    using GridDiscretization = GetPropType<SubDomainTypeTag<id>, Properties::GridGeometry>;

    template<std::size_t id>
    using Problem = typename MDTraits::template SubDomain<id>::Problem;

    using JacobianMatrix = typename MDTraits::JacobianMatrix;
    using SolutionVector = typename MDTraits::SolutionVector;
    using ResidualType = typename MDTraits::ResidualVector;
    using CouplingManager = CMType;

private:
    using ProblemTuple = typename MDTraits::template TupleOfSharedPtrConst<Problem>;
    using GridDiscretizationTuple = typename MDTraits::template TupleOfSharedPtrConst<GridDiscretization>;
    using GridVariablesTuple = typename MDTraits::template TupleOfSharedPtr<GridVariables>;
    using ThisType = MultiStageMultiDomainAssembler<MDTraits, CouplingManager, diffMethod>;

    template<std::size_t id>
    using SubDomainAssemblerView = MultiDomainAssemblerSubDomainView<ThisType, id>;

    template<class DiscretizationMethod, std::size_t id>
    struct SubDomainAssemblerType;

    template<std::size_t id, class DM>
    struct SubDomainAssemblerType<DiscretizationMethods::CVFE<DM>, id>
    {
        using type = Dumux::Experimental::SubDomainCVFELocalAssembler<id, SubDomainTypeTag<id>, SubDomainAssemblerView<id>, diffMethod, true>;
    };

    template<std::size_t id>
    using SubDomainAssembler = typename SubDomainAssemblerType<typename GridDiscretization<id>::DiscretizationMethod, id>::type;

public:
    /*!
     * \brief The constructor
     *
     * \param problem The problems of all subdomains
     * \param gridDiscretization The grid discretizations of all subdomains
     * \param gridVariables The grid variables of all subdomains
     * \param couplingManager The coupling manager
     * \param msMethod The multi-stage time stepping method
     * \param prevSol The solution at the previous time level
     */
    MultiStageMultiDomainAssembler(ProblemTuple problem,
                                   GridDiscretizationTuple gridDiscretization,
                                   GridVariablesTuple gridVariables,
                                   std::shared_ptr<CouplingManager> couplingManager,
                                   std::shared_ptr<const Experimental::MultiStageMethod<Scalar>> msMethod,
                                   const SolutionVector& prevSol)
    : couplingManager_(couplingManager)
    , timeSteppingMethod_(msMethod)
    , problemTuple_(std::move(problem))
    , gridDiscretizationTuple_(std::move(gridDiscretization))
    , gridVariablesTuple_(std::move(gridVariables))
    , prevSol_(&prevSol)
    {
        std::cout << "Instantiated multi-stage multi-domain assembler." << std::endl;

        enableMultithreading_ = CouplingManagerSupportsMultithreadedAssembly<CouplingManager>::value
            && Grid::Capabilities::allGridsSupportsMultithreading(gridDiscretizationTuple_)
            && !Multithreading::isSerial()
            && getParam<bool>("Assembly.Multithreading", true);

        maybeComputeColors_();
    }

    /*!
     * \brief Assembles the Jacobian and residual of the current stage
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

            // add the contributions of all previous stages
            auto constantResidualComponent = (*residual_)[domainId];
            constantResidualComponent = 0.0;
            for (std::size_t k = 0; k < stageParams_->size()-1; ++k)
            {
                if (!stageParams_->skipTemporal(k))
                    constantResidualComponent.axpy(stageParams_->temporalWeight(k), temporalOperatorEvaluations_[k][domainId]);
                if (!stageParams_->skipSpatial(k))
                    constantResidualComponent.axpy(stageParams_->spatialWeight(k), spatialOperatorEvaluations_[k][domainId]);
            }

            for (std::size_t i = 0; i < constantResidualComponent.size(); ++i)
                for (std::size_t ii = 0; ii < constantResidualComponent[i].size(); ++ii)
                    (*residual_)[domainId][i][ii] += constrainedDofs_[domainId][i][ii] > 0.5 ? 0.0 : constantResidualComponent[i][ii];

            enforceProblemConstraints_(domainId, (*jacobian_)[domainId], (*residual_)[domainId], curSol[domainId]);
        });
    }

    //! Assembling the residual without the Jacobian is not supported
    void assembleResidual(const SolutionVector&)
    { DUNE_THROW(Dune::NotImplemented, "assembleResidual without Jacobian for multi-stage multi-domain assembly"); }

    //! Set up the Jacobian and residual storage
    void setLinearSystem()
    {
        jacobian_ = std::make_shared<JacobianMatrix>();
        residual_ = std::make_shared<ResidualType>();

        setJacobianBuildMode_(*jacobian_);
        setJacobianPattern_(*jacobian_);
        setResidualSize_(*residual_);
    }

    //! Update the grid variables and the coupling manager for the given solution
    void updateGridVariables(const SolutionVector& curSol)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(gridVariablesTuple_)), [&](const auto domainId)
        { this->gridVariables(domainId).update(curSol[domainId]); });

        // cross-domain quantities have to reflect the current iterate
        couplingManager_->updateSolution(curSol);
    }

    //! Reset the grid variables to the given solution (e.g. after a failed time step)
    void resetTimeStep(const SolutionVector& curSol)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(gridVariablesTuple_)), [&](const auto domainId)
        { this->gridVariables(domainId).resetTimeStep(curSol[domainId]); });

        this->clearStages();
    }

    //! Prepare the assembler for the next stage of the time integration
    template<class StageParamsPtr>
    void prepareStage(SolutionVector& x, StageParamsPtr params)
    {
        stageParams_ = std::move(params);
        const auto curStage = stageParams_->size() - 1;

        // in the first stage, also assemble the operators at the previous time level (stage 0)
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
            assembleOperatorEvaluations_(x);
        }

        // The evaluations recorded while solving the previous stage belong to the last iterate
        // the solver assembled, which is not the stage solution for a solver that stops after an
        // update (e.g. Newton on the shift criterion) or assembles only once (a linear solver).
        else
        {
            if (spatialOperatorEvaluations_.size() != curStage)
                DUNE_THROW(Dune::InvalidStateException, "Invalid state. Maybe you forgot to call clearStages()");

            using namespace Dune::Hybrid;
            forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
            {
                setProblemTime_(*std::get<domainId>(problemTuple_), stageParams_->timeAtStage(curStage-1));
            });
            assembleOperatorEvaluations_(x);
        }

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
        {
            setProblemTime_(*std::get<domainId>(problemTuple_), stageParams_->timeAtStage(curStage));
        });

        resetResidual_();
        spatialOperatorEvaluations_.push_back(*residual_);
        temporalOperatorEvaluations_.push_back(*residual_);
    }

    //! Clear all stored stage evaluations
    void clearStages()
    {
        spatialOperatorEvaluations_.clear();
        temporalOperatorEvaluations_.clear();
        stageParams_.reset();
    }

    bool isStationaryProblem() const
    { return false; }

    bool isImplicit() const
    { return timeSteppingMethod_->implicit(); }

    //! the number of dof locations of domain i
    template<std::size_t i>
    std::size_t numDofs(Dune::index_constant<i> domainId) const
    { return std::get<domainId>(gridDiscretizationTuple_)->numDofs(); }

    //! the problem of domain i
    template<std::size_t i>
    const auto& problem(Dune::index_constant<i> domainId) const
    { return *std::get<domainId>(problemTuple_); }

    //! the grid discretization of domain i
    template<std::size_t i>
    const auto& gridDiscretization(Dune::index_constant<i> domainId) const
    { return *std::get<domainId>(gridDiscretizationTuple_); }

    //! the grid view of domain i
    template<std::size_t i>
    const auto& gridView(Dune::index_constant<i> domainId) const
    { return gridDiscretization(domainId).gridView(); }

    //! the grid variables of domain i
    template<std::size_t i>
    GridVariables<i>& gridVariables(Dune::index_constant<i> domainId)
    { return *std::get<domainId>(gridVariablesTuple_); }

    //! the grid variables of domain i
    template<std::size_t i>
    const GridVariables<i>& gridVariables(Dune::index_constant<i> domainId) const
    { return *std::get<domainId>(gridVariablesTuple_); }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    JacobianMatrix& jacobian()
    { return *jacobian_; }

    ResidualType& residual()
    { return *residual_; }

    const SolutionVector& prevSol() const
    { return *prevSol_; }

    //! Set the solution at the previous time level
    void setPreviousSolution(const SolutionVector& sol)
    { prevSol_ = &sol; }

    //! The temporal and spatial weight of the current stage
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
    //! Assemble the unweighted temporal and spatial operators at x into the last stored evaluations
    void assembleOperatorEvaluations_(const SolutionVector& x)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
        {
            auto& spatial = spatialOperatorEvaluations_.back()[domainId];
            auto& temporal = temporalOperatorEvaluations_.back()[domainId];
            spatial = 0.0;
            temporal = 0.0;
            assemble_(domainId, [&](const auto& element)
            {
                MultiDomainAssemblerSubDomainView view{*this, domainId};
                SubDomainAssembler<domainId()> subDomainAssembler(view, element, x, *couplingManager_);
                subDomainAssembler.assembleCurrentResidual(temporal, spatial);
            });
        });
    }
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

    //! Enforce the Dirichlet constraints of problem.constraints() in the Jacobian row and residual of domain i
    template<std::size_t i, class JacRow, class Res, class Sol>
    void enforceProblemConstraints_(Dune::index_constant<i> domainI, JacRow& jacRow, Res& res, const Sol& curSol) const
    {
        if constexpr (Dumux::Detail::hasSubProblemGlobalConstraints<Problem<domainI>>())
        {
            auto& jac = jacRow[domainI];

            auto applyDirichletConstraint = [&](const auto& dofIdx, const auto& values,
                                                const auto eqIdx, const auto pvIdx)
            {
                res[dofIdx][eqIdx] = curSol[dofIdx][pvIdx] - values[pvIdx];

                auto& row = jac[dofIdx];
                for (auto col = row.begin(); col != row.end(); ++col)
                    row[col.index()][eqIdx] = 0.0;
                jac[dofIdx][dofIdx][eqIdx][pvIdx] = 1.0;

                using namespace Dune::Hybrid;
                forEach(makeIncompleteIntegerSequence<JacobianMatrix::N(), domainI>(), [&](const auto couplingId)
                {
                    auto& rowCoupling = jacRow[couplingId][dofIdx];
                    for (auto c = rowCoupling.begin(); c != rowCoupling.end(); ++c)
                        rowCoupling[c.index()][eqIdx] = 0.0;
                });
            };

            for (const auto& constraintData : this->problem(domainI).constraints())
            {
                const auto& info = constraintData.constraintInfo();
                const auto& values = constraintData.values();
                const auto dofIdx = constraintData.dofIndex();
                for (int eqIdx = 0; eqIdx < info.size(); ++eqIdx)
                    if (info.isConstraintEquation(eqIdx))
                        applyDirichletConstraint(dofIdx, values, eqIdx, info.eqToPriVarIndex(eqIdx));
            }
        }
    }

    template<std::size_t i, std::size_t j, typename std::enable_if_t<(i==j), int> = 0>
    Dune::MatrixIndexSet getJacobianPattern_(Dune::index_constant<i> domainI,
                                             Dune::index_constant<j>) const
    {
        const auto& gridDisc = gridDiscretization(domainI);
        auto pattern = timeSteppingMethod_->implicit() ? getJacobianPattern<true>(gridDisc)
                                                       : getJacobianPattern<false>(gridDisc);
        couplingManager_->extendJacobianPattern(domainI, pattern);
        return pattern;
    }

    template<std::size_t i, std::size_t j, typename std::enable_if_t<(i!=j), int> = 0>
    Dune::MatrixIndexSet getJacobianPattern_(Dune::index_constant<i> domainI,
                                             Dune::index_constant<j> domainJ) const
    {
        if (timeSteppingMethod_->implicit())
            return getCouplingJacobianPattern<true>(*couplingManager_,
                domainI, gridDiscretization(domainI),
                domainJ, gridDiscretization(domainJ));
        else
            return getCouplingJacobianPattern<false>(*couplingManager_,
                domainI, gridDiscretization(domainI),
                domainJ, gridDiscretization(domainJ));
    }

    template<class P>
    void setProblemTime_(const P& p, const Scalar t)
    {
        if constexpr (requires { p.setTime(t); })
            p.setTime(t);
        else
            static_assert(!requires (P& q) { q.setTime(Scalar{}); },
                "The multi-stage assembler sets the stage time through a const problem: setTime has to be const.");
    }

    std::shared_ptr<const Experimental::MultiStageMethod<Scalar>> timeSteppingMethod_;
    std::vector<ResidualType> spatialOperatorEvaluations_;
    std::vector<ResidualType> temporalOperatorEvaluations_;
    ResidualType constrainedDofs_;
    std::shared_ptr<const StageParams> stageParams_;

    ProblemTuple problemTuple_;
    GridDiscretizationTuple gridDiscretizationTuple_;
    GridVariablesTuple gridVariablesTuple_;

    const SolutionVector* prevSol_ = nullptr;

    std::shared_ptr<JacobianMatrix> jacobian_;
    std::shared_ptr<ResidualType> residual_;

    bool enableMultithreading_ = false;
};

} // end namespace Dumux::Experimental

#endif
