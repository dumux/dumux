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
 * \brief A linear system assembler for the implicit solves of time integration methods
 *        for multiple domains discretized with the control-volume finite element interface.
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
#include <optional>
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

#include <dumux/experimental/timestepping/timelevel.hh>
#include <dumux/experimental/timestepping/timeintegrationmethod.hh>
#include <dumux/experimental/timestepping/implicitsolve.hh>

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
    using SolveParams = Experimental::ImplicitSolveParams<Scalar, ResidualType>;
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
     * \param method The time integration method, which determines the Jacobian pattern
     * \param prevSol The solution at the previous time level (used by solvers that reset a time step)
     */
    MultiStageMultiDomainAssembler(ProblemTuple problem,
                                   GridDiscretizationTuple gridDiscretization,
                                   GridVariablesTuple gridVariables,
                                   std::shared_ptr<CouplingManager> couplingManager,
                                   std::shared_ptr<const Experimental::TimeIntegrationMethod<Scalar>> method,
                                   const SolutionVector& prevSol)
    : couplingManager_(couplingManager)
    , timeSteppingMethod_(method)
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
     * \brief Assembles the Jacobian and residual of the current implicit solve
     */
    void assembleJacobianAndResidual(const SolutionVector& curSol)
    {
        if (!solveParams_)
            DUNE_THROW(Dune::InvalidStateException, "No implicit solve set. Call setImplicitSolve first.");

        resetJacobian_();
        resetResidual_();

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
        {
            auto& jacRow = (*jacobian_)[domainId];
            assemble_(domainId, [&](const auto& element)
            {
                MultiDomainAssemblerSubDomainView view{*this, domainId};
                SubDomainAssembler<domainId()> subDomainAssembler(view, element, curSol, *couplingManager_);
                subDomainAssembler.assembleJacobianAndResidual(
                    jacRow, (*residual_)[domainId],
                    gridVariablesTuple_,
                    solveParams_->weights,
                    constrainedDofs_[domainId]
                );
            });

            // the contribution of known states does not enter the constraint equations
            if (solveParams_->constant)
            {
                const auto& constant = (*solveParams_->constant)[domainId];
                for (std::size_t i = 0; i < constant.size(); ++i)
                    for (std::size_t ii = 0; ii < constant[i].size(); ++ii)
                        (*residual_)[domainId][i][ii] += constrainedDofs_[domainId][i][ii] > 0.5 ? 0.0 : constant[i][ii];
            }

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
    }

    /*!
     * \brief Set the weights, the time level and the contribution of known states for the next implicit solve
     */
    void setImplicitSolve(const SolutionVector&, const SolveParams& params)
    {
        if (params.weights.inertia != 0.0)
            DUNE_THROW(Dune::NotImplemented, "Inertia terms in the multi-stage multi-domain assembler");

        setProblemTimes_(params.timeLevel.current());
        solveParams_ = params;
    }

    /*!
     * \brief Evaluate the unweighted operators at the given solution and time level
     * \note The grid variables and the coupling manager have to be updated to the given solution.
     */
    void assembleOperators(const SolutionVector& sol,
                           const Experimental::TimeLevel<Scalar>& timeLevel,
                           Experimental::OperatorEvaluations<ResidualType> evaluations)
    {
        if (evaluations.inertia)
            DUNE_THROW(Dune::NotImplemented, "Inertia terms in the multi-stage multi-domain assembler");

        setProblemTimes_(timeLevel.current());

        ResidualType unused;
        auto& storage = evaluations.storage ? *evaluations.storage : unused;
        auto& spatial = evaluations.spatial ? *evaluations.spatial : unused;
        for (auto* r : {&storage, &spatial})
        {
            setResidualSize_(*r);
            (*r) = 0.0;
        }

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
        {
            assemble_(domainId, [&](const auto& element)
            {
                MultiDomainAssemblerSubDomainView view{*this, domainId};
                SubDomainAssembler<domainId()> subDomainAssembler(view, element, sol, *couplingManager_);
                subDomainAssembler.assembleCurrentResidual(storage[domainId], spatial[domainId]);
            });
        });
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

    //! The operator weights of the current implicit solve
    const Experimental::OperatorWeights<Scalar>& currentOperatorWeights() const
    {
        if (!solveParams_)
            DUNE_THROW(Dune::InvalidStateException, "No implicit solve set. Call setImplicitSolve first.");
        return solveParams_->weights;
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
    void setProblemTimes_(const Scalar t)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
        { setProblemTime_(*std::get<domainId>(problemTuple_), t); });
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

    std::shared_ptr<const Experimental::TimeIntegrationMethod<Scalar>> timeSteppingMethod_;
    ResidualType constrainedDofs_;
    std::optional<SolveParams> solveParams_;

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
