// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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

#include <dumux/experimental/assembly/subdomaincclocalassembler.hh>
#include <dumux/experimental/assembly/subdomaincvfelocalassembler.hh>
#include <dumux/experimental/assembly/multistagefvlocaloperator.hh>

#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>

namespace Dumux::Grid::Capabilities {

namespace Detail {
// helper for multi-domain models
template<class T, std::size_t... I>
bool allGridsSupportsMultithreadingImpl(const T& gridGeometries, std::index_sequence<I...>)
{
    return (... && supportsMultithreading(std::get<I>(gridGeometries)->gridView()));
}
} // end namespace Detail

// helper for multi-domain models (all grids have to support multithreading)
template<class... GG>
bool allGridsSupportsMultithreading(const std::tuple<GG...>& gridGeometries)
{
    return Detail::allGridsSupportsMultithreadingImpl<std::tuple<GG...>>(gridGeometries, std::make_index_sequence<sizeof...(GG)>());
}

} // end namespace Dumux::Grid::Capabilities

namespace Dumux {

/*!
 * \ingroup Experimental
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief Type trait that is specialized for coupling manager supporting multithreaded assembly
 * \note A coupling manager implementation that wants to enable multithreaded assembly has to specialize this trait
 */
template<class CM>
struct CouplingManagerSupportsMultithreadedAssembly : public std::false_type
{};

} // end namespace Dumux

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
    using StageParams = Experimental::MultiStageParams<Scalar>;

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
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    MultiStageMultiDomainFVAssembler(ProblemTuple problem,
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
        resetJacobian_();

        resetResidual_();
        spatialOperatorEvaluations_.back() = 0.0;
        temporalOperatorEvaluations_.back() = 0.0;

        if (stageParams_->size() != spatialOperatorEvaluations_.size())
            DUNE_THROW(Dune::InvalidStateException, "Wrong number of residuals");

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
        {
            // assemble the spatial and temporal residual of the current time step and the Jacobian
            // w.r.t to the current solution (the current solution on the current stage)
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

            // assemble the full residual for the time integration stage
            auto constantResidualComponent = (*residual_)[domainId];
            constantResidualComponent = 0.0;
            for (std::size_t k = 0; k < stageParams_->size()-1; ++k)
            {
                if (!stageParams_->skipTemporal(k))
                    constantResidualComponent.axpy(stageParams_->temporalWeight(k), temporalOperatorEvaluations_[k][domainId]);
                if (!stageParams_->skipSpatial(k))
                    constantResidualComponent.axpy(stageParams_->spatialWeight(k), spatialOperatorEvaluations_[k][domainId]);
            }

            // masked summation of constant residual component onto this stage's resiudal component
            for (std::size_t i = 0; i < constantResidualComponent.size(); ++i)
                for (std::size_t ii = 0; ii < constantResidualComponent[i].size(); ++ii)
                    (*residual_)[domainId][i][ii] += constrainedDofs_[domainId][i][ii] > 0.5 ? 0.0 : constantResidualComponent[i][ii];
        });
    }

    //! compute the residuals using the internal residual
    void assembleResidual(const SolutionVector& curSol)
    { DUNE_THROW(Dune::NotImplemented, "residual"); }

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

    void clearStages()
    {
        spatialOperatorEvaluations_.clear();
        temporalOperatorEvaluations_.clear();
        stageParams_.reset();
    }

    template<class StageParams>
    void prepareStage(SolutionVector& x, StageParams params)
    {
        stageParams_ = std::move(params);
        const auto curStage = stageParams_->size() - 1;

        // in the first stage, also assemble the old residual
        if (curStage == 1)
        {
            // update time in variables?
            using namespace Dune::Hybrid;
            forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
            {
                setProblemTime_(*std::get<domainId>(problemTuple_), stageParams_->timeAtStage(curStage));
            });

            resetResidual_(); // residual resized and zero
            spatialOperatorEvaluations_.push_back(*residual_);
            temporalOperatorEvaluations_.push_back(*residual_);

            // assemble stage 0 residuals
            using namespace Dune::Hybrid;
            forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
            {
                auto& spatial = spatialOperatorEvaluations_.back()[domainId];
                auto& temporal = temporalOperatorEvaluations_.back()[domainId];
                assemble_(domainId, [&](const auto& element)
                {
                    MultiDomainAssemblerSubDomainView view{*this, domainId};
                    SubDomainAssembler<domainId()> subDomainAssembler(view, element, x, *couplingManager_);
                    subDomainAssembler.localResidual().spatialWeight(1.0);
                    subDomainAssembler.localResidual().temporalWeight(1.0);
                    subDomainAssembler.assembleCurrentResidual(spatial, temporal);
                });
            });
        }

        // update time in variables?
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
        {
            setProblemTime_(*std::get<domainId>(problemTuple_), stageParams_->timeAtStage(curStage));
        });

        resetResidual_(); // residual resized and zero
        spatialOperatorEvaluations_.push_back(*residual_);
        temporalOperatorEvaluations_.push_back(*residual_);
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

    // TODO make this nicer with a is_detected trait in a common location
    template<class P>
    void setProblemTime_(const P& p, const Scalar t)
    { setProblemTimeImpl_(p, t, 0); }

    template<class P>
    auto setProblemTimeImpl_(const P& p, const Scalar t, int) -> decltype(p.setTime(0))
    { p.setTime(t); }

    template<class P>
    void setProblemTimeImpl_(const P& p, const Scalar t, long)
    {}

    std::shared_ptr<const Experimental::MultiStageMethod<Scalar>> timeSteppingMethod_;
    std::vector<ResidualType> spatialOperatorEvaluations_;
    std::vector<ResidualType> temporalOperatorEvaluations_;
    ResidualType constrainedDofs_;
    std::shared_ptr<const StageParams> stageParams_;

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
