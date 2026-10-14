// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for finite volume schemes
 *        with multiple domains
 */
#ifndef DUMUX_EXPERIMENTAL_MULTISTAGE_MULTIDOMAIN_FV_ASSEMBLER_HH
#define DUMUX_EXPERIMENTAL_MULTISTAGE_MULTIDOMAIN_FV_ASSEMBLER_HH

#include <vector>
#include <type_traits>
#include <tuple>
#include <memory>
#include <optional>

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
#include <dumux/multidomain/assemblytraits.hh>

#include <dumux/experimental/assembly/subdomaincclocalassembler.hh>
#include <dumux/experimental/assembly/subdomaincvfelocalassembler.hh>
#include <dumux/experimental/assembly/multistagefvlocaloperator.hh>

#include <dumux/experimental/timestepping/timelevel.hh>
#include <dumux/experimental/timestepping/timeintegrationmethod.hh>
#include <dumux/experimental/timestepping/implicitsolve.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for finite volume schemes (box, tpfa, mpfa, ...)
 *        with multiple domains
 * \tparam MDTraits the multidimensional traits
 * \tparam diffMethod the differentiation method to residual compute derivatives
 */
template<class MDTraits, class CMType, DiffMethod diffMethod>
class MultiStageMultiDomainFVAssembler
{
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

public:
    using Traits = MDTraits;

    using Scalar = typename MDTraits::Scalar;

    //! TODO get rid of this GetPropType
    template<std::size_t id>
    using LocalResidual = GetPropType<SubDomainTypeTag<id>, Properties::LocalResidual>;

    template<std::size_t id>
    using GridVariables = typename MDTraits::template SubDomain<id>::GridVariables;

    template<std::size_t id>
    using GridGeometry = typename MDTraits::template SubDomain<id>::GridGeometry;

    template<std::size_t id>
    using Problem = typename MDTraits::template SubDomain<id>::Problem;

    using JacobianMatrix = typename MDTraits::JacobianMatrix;
    using SolutionVector = typename MDTraits::SolutionVector;
    using ResidualType = typename MDTraits::ResidualVector;
    using SolveParams = Experimental::ImplicitSolveParams<Scalar, ResidualType>;

    using CouplingManager = CMType;

private:

    using ProblemTuple = typename MDTraits::template TupleOfSharedPtrConst<Problem>;
    using GridGeometryTuple = typename MDTraits::template TupleOfSharedPtrConst<GridGeometry>;
    using GridVariablesTuple = typename MDTraits::template TupleOfSharedPtr<GridVariables>;

    using ThisType = MultiStageMultiDomainFVAssembler<MDTraits, CouplingManager, diffMethod>;

    template<std::size_t id>
    using SubDomainAssemblerView = MultiDomainAssemblerSubDomainView<ThisType, id>;

    template<class DiscretizationMethod, std::size_t id>
    struct SubDomainAssemblerType;

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethods::CCTpfa, id>
    {
        using type = Experimental::SubDomainCCLocalAssembler<id, SubDomainTypeTag<id>, SubDomainAssemblerView<id>, diffMethod>;
    };

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethods::CCMpfa, id>
    {
        using type = Experimental::SubDomainCCLocalAssembler<id, SubDomainTypeTag<id>, SubDomainAssemblerView<id>, diffMethod>;
    };

    template<std::size_t id, class DM>
    struct SubDomainAssemblerType<DiscretizationMethods::CVFE<DM>, id>
    {
        using type = Experimental::SubDomainCVFELocalAssembler<id, SubDomainTypeTag<id>, SubDomainAssemblerView<id>, diffMethod>;
    };

    template<std::size_t id>
    using SubDomainAssembler = typename SubDomainAssemblerType<typename GridGeometry<id>::DiscretizationMethod, id>::type;

public:
    /*!
     * \brief The constructor for instationary problems
     * \param method The time integration method, which determines the Jacobian pattern
     * \param prevSol The solution at the beginning of the time step (used by solvers that reset a time step)
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    MultiStageMultiDomainFVAssembler(ProblemTuple problem,
                                     GridGeometryTuple gridGeometry,
                                     GridVariablesTuple gridVariables,
                                     std::shared_ptr<CouplingManager> couplingManager,
                                     std::shared_ptr<const Experimental::TimeIntegrationMethod<Scalar>> method,
                                     const SolutionVector& prevSol)
    : couplingManager_(couplingManager)
    , timeSteppingMethod_(method)
    , problemTuple_(std::move(problem))
    , gridGeometryTuple_(std::move(gridGeometry))
    , gridVariablesTuple_(std::move(gridVariables))
    , prevSol_(&prevSol)
    {
        std::cout << "Instantiated assembler for an instationary problem." << std::endl;

        enableMultithreading_ = CouplingManagerSupportsMultithreadedAssembly<CouplingManager>::value
            && Grid::Capabilities::allGridsSupportsMultithreading(gridGeometryTuple_)
            && !Multithreading::isSerial()
            && getParam<bool>("Assembly.Multithreading", true);

        maybeComputeColors_();
    }

    /*!
     * \brief Assembles the global Jacobian of the residual
     *        and the residual for the current solution.
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
        });
    }

    /*!
     * \brief Assembles the residual for the given solution without the Jacobian.
     *
     * Required by everything in the Newton solver that evaluates the residual somewhere
     * other than the linearization point: a line search, and a residual-based convergence
     * criterion. Without it neither can be used with a multi-stage method.
     *
     * The residual is the weighted sum of the storage and spatial operators at `curSol`
     * plus the contribution of known states set with setImplicitSolve.
     *
     * The Dirichlet constraints of control-volume finite element subdomains are evaluated here
     * rather than taken from the mask of a previous Jacobian assembly, which need not have
     * happened in this stage. The residual of a constrained degree of freedom is the constraint
     * itself and does not depend on the stage weights, so it overwrites the weighted sum.
     */
    void assembleResidual(const SolutionVector& curSol)
    {
        if (!solveParams_)
            DUNE_THROW(Dune::InvalidStateException, "No implicit solve set. Call setImplicitSolve first.");

        if (!residual_)
        {
            residual_ = std::make_shared<ResidualType>();
            setResidualSize_(*residual_);
            setResidualSize_(constrainedDofs_);
            constrainedDofs_ = 0.0;
        }

        const auto& weights = solveParams_->weights;
        ResidualType storage, spatial;
        assembleOperators_(curSol, {.storage = &storage, .spatial = &spatial});

        (*residual_) = 0.0;
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
        {
            auto& residual = (*residual_)[domainId];
            if (weights.storage != 0.0)
                residual.axpy(weights.storage, storage[domainId]);
            if (weights.spatial != 0.0)
                residual.axpy(weights.spatial, spatial[domainId]);
            if (solveParams_->constant)
                residual += (*solveParams_->constant)[domainId];

            constrainedDofs_[domainId] = 0.0;
            if constexpr (DiscretizationMethods::isCVFE<typename GridGeometry<domainId>::DiscretizationMethod>)
            {
                assemble_(domainId, [&](const auto& element)
                {
                    MultiDomainAssemblerSubDomainView view{*this, domainId};
                    SubDomainAssembler<domainId()> subDomainAssembler(view, element, curSol, *couplingManager_);
                    subDomainAssembler.bindLocalViews();
                    subDomainAssembler.enforceDirichletConstraints(
                        [&](const auto& scvI, const auto& dirichletValues, const auto eqIdx, const auto pvIdx)
                        {
                            residual[scvI.dofIndex()][eqIdx]
                                = subDomainAssembler.curElemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];
                            constrainedDofs_[domainId][scvI.dofIndex()][eqIdx] = 1.0;
                        });
                });
            }
        });
    }

    /*!
     * \brief The version without arguments uses the default constructor to create
     *        the jacobian and residual objects in this assembler if you don't need them outside this class
     */
    void setLinearSystem()
    {
        jacobian_ = std::make_shared<JacobianMatrix>();
        residual_ = std::make_shared<ResidualType>();

        setJacobianBuildMode_(*jacobian_);
        setJacobianPattern_(*jacobian_);
        setResidualSize_(*residual_);
    }

    /*!
     * \brief Updates the grid variables with the given solution
     */
    void updateGridVariables(const SolutionVector& curSol)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(gridVariablesTuple_)), [&](const auto domainId)
        { this->gridVariables(domainId).update(curSol[domainId]); });
    }

    /*!
     * \brief Resets the grid variables to the last time step
     */
    void resetTimeStep(const SolutionVector& curSol)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(gridVariablesTuple_)), [&](const auto domainId)
        { this->gridVariables(domainId).resetTimeStep(curSol[domainId]); });
    }

    //! the number of dof locations of domain i
    template<std::size_t i>
    std::size_t numDofs(Dune::index_constant<i> domainId) const
    { return std::get<domainId>(gridGeometryTuple_)->numDofs(); }

    //! the problem of domain i
    template<std::size_t i>
    const auto& problem(Dune::index_constant<i> domainId) const
    { return *std::get<domainId>(problemTuple_); }

    //! the finite volume grid geometry of domain i
    template<std::size_t i>
    const auto& gridGeometry(Dune::index_constant<i> domainId) const
    { return *std::get<domainId>(gridGeometryTuple_); }

    //! the grid view of domain i
    template<std::size_t i>
    const auto& gridView(Dune::index_constant<i> domainId) const
    { return gridGeometry(domainId).gridView(); }

    //! the grid variables of domain i
    template<std::size_t i>
    GridVariables<i>& gridVariables(Dune::index_constant<i> domainId)
    { return *std::get<domainId>(gridVariablesTuple_); }

    //! the grid variables of domain i
    template<std::size_t i>
    const GridVariables<i>& gridVariables(Dune::index_constant<i> domainId) const
    { return *std::get<domainId>(gridVariablesTuple_); }

    //! the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    //! the full Jacobian matrix
    JacobianMatrix& jacobian()
    { return *jacobian_; }

    //! the full residual vector
    ResidualType& residual()
    { return *residual_; }

    //! the solution before time integration
    const SolutionVector& prevSol() const
    { return *prevSol_; }

    /*!
     * \brief Create a local residual object (used by the local assembler)
     */
    template<std::size_t i>
    MultiStageFVLocalOperator<LocalResidual<i>> localResidual(Dune::index_constant<i> domainId) const
    { return { LocalResidual<i>{std::get<domainId>(problemTuple_).get(), nullptr} }; }

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
        assembleOperators_(sol, evaluations);
    }

    //! TODO get rid of this (called by Newton but shouldn't be necessary)
    bool isStationaryProblem() const
    { return false; }

    bool isImplicit() const
    { return timeSteppingMethod_->implicit(); }

protected:
    //! the coupling manager coupling the sub domains
    std::shared_ptr<CouplingManager> couplingManager_;

private:
    /*!
     * \brief Sets the jacobian build mode
     */
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
                    DUNE_THROW(Dune::NotImplemented, "Only BCRS matrices with random build mode are supported at the moment");
            });
        });
    }

    /*!
     * \brief Sets the jacobian sparsity pattern.
     */
    void setJacobianPattern_(JacobianMatrix& jac) const
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainI)
        {
            forEach(integralRange(Dune::Hybrid::size(jac[domainI])), [&](const auto domainJ)
            {
                const auto pattern = this->getJacobianPattern_(domainI, domainJ);
                pattern.exportIdx(jac[domainI][domainJ]);
            });
        });
    }

    /*!
     * \brief Resizes the residual
     */
    void setResidualSize_(ResidualType& res) const
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(res)), [&](const auto domainId)
        { res[domainId].resize(this->numDofs(domainId)); });
    }

    // reset the residual vector to 0.0
    void resetResidual_()
    {
        if(!residual_)
        {
            residual_ = std::make_shared<ResidualType>();
            setResidualSize_(*residual_);
        }

        setResidualSize_(constrainedDofs_);

        (*residual_) = 0.0;
        constrainedDofs_ = 0.0;
    }

    // reset the jacobian vector to 0.0
    void resetJacobian_()
    {
        if(!jacobian_)
        {
            jacobian_ = std::make_shared<JacobianMatrix>();
            setJacobianBuildMode_(*jacobian_);
            setJacobianPattern_(*jacobian_);
        }

       (*jacobian_)  = 0.0;
    }

    //! Computes the colors
    void maybeComputeColors_()
    {
        if constexpr (CouplingManagerSupportsMultithreadedAssembly<CouplingManager>::value)
            if (enableMultithreading_)
                couplingManager_->computeColorsForAssembly();
    }

    template<std::size_t i, class SubRes>
    void assembleResidual_(Dune::index_constant<i> domainId, SubRes& subRes,
                           const SolutionVector& curSol)
    {
        DUNE_THROW(Dune::NotImplemented, "assembleResidual_");
    }

    /*!
     * \brief A method assembling something per element
     * \note Handles exceptions for parallel runs
     * \throws NumericalProblem on all processes if an exception is thrown during assembly
     */
    template<std::size_t i, class AssembleElementFunc>
    void assemble_(Dune::index_constant<i> domainId, AssembleElementFunc&& assembleElement) const
    {
        // a state that will be checked on all processes
        bool succeeded = false;

        // try assembling using the local assembly function
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

            // fallback for coupling managers that don't support multithreaded assembly (yet)
            // or if multithreaded assembly is disabled
            // let the local assembler add the element contributions
            for (const auto& element : elements(gridView(domainId)))
                assembleElement(element);

            // if we get here, everything worked well on this process
            succeeded = true;
        }
        // throw exception if a problem occurred
        catch (NumericalProblem &e)
        {
            std::cout << "rank " << gridView(domainId).comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n";
            succeeded = false;
        }

        // make sure everything worked well on all processes
        if (gridView(domainId).comm().size() > 1)
            succeeded = gridView(domainId).comm().min(succeeded);

        // if not succeeded rethrow the error on all processes
        if (!succeeded)
            DUNE_THROW(NumericalProblem, "A process did not succeed in linearizing the system");
    }

    // get diagonal block pattern
    template<std::size_t i, std::size_t j, typename std::enable_if_t<(i==j), int> = 0>
    Dune::MatrixIndexSet getJacobianPattern_(Dune::index_constant<i> domainI,
                                             Dune::index_constant<j> domainJ) const
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

    // get coupling block pattern
    template<std::size_t i, std::size_t j, typename std::enable_if_t<(i!=j), int> = 0>
    Dune::MatrixIndexSet getJacobianPattern_(Dune::index_constant<i> domainI,
                                             Dune::index_constant<j> domainJ) const
    {
        if (timeSteppingMethod_->implicit())
            return getCouplingJacobianPattern<true>(*couplingManager_,
                domainI, gridGeometry(domainI),
                domainJ, gridGeometry(domainJ)
            );
        else
            return getCouplingJacobianPattern<false>(*couplingManager_,
                domainI, gridGeometry(domainI),
                domainJ, gridGeometry(domainJ)
            );
    }

    //! Evaluate the unweighted operators at the given solution at the current problem time
    void assembleOperators_(const SolutionVector& sol, Experimental::OperatorEvaluations<ResidualType> evaluations)
    {
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
                subDomainAssembler.localResidual().spatialWeight(1.0, !evaluations.spatial);
                subDomainAssembler.localResidual().temporalWeight(1.0, !evaluations.storage);
                subDomainAssembler.assembleCurrentResidual(spatial[domainId], storage[domainId]);
            });
        });
    }

    void setProblemTimes_(const Scalar t)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
        { setProblemTime_(*std::get<domainId>(problemTuple_), t); });
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

    //! pointer to the problem to be solved
    ProblemTuple problemTuple_;

    //! the finite volume geometry of the grid
    GridGeometryTuple gridGeometryTuple_;

    //! the variables container for the grid
    GridVariablesTuple gridVariablesTuple_;

    //! Pointer to the previous solution
    const SolutionVector* prevSol_;

    //! shared pointers to the jacobian matrix and residual
    std::shared_ptr<JacobianMatrix> jacobian_;
    std::shared_ptr<ResidualType> residual_;

    //! if multithreaded assembly is enabled
    bool enableMultithreading_ = false;
};

} // end namespace Dumux

#endif
