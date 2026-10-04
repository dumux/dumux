// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \brief A multi-stage (Runge-Kutta) linear system assembler for general discretization
 *        schemes, built on the localDof/grid-variables local assembler used by Assembler.
 */
#ifndef DUMUX_EXPERIMENTAL_MULTISTAGE_ASSEMBLER_HH
#define DUMUX_EXPERIMENTAL_MULTISTAGE_ASSEMBLER_HH

#include <vector>
#include <deque>
#include <memory>
#include <utility>

#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/gridcapabilities.hh>
#include <dumux/common/typetraits/vector.hh>

#include <dumux/discretization/method.hh>
#include <dumux/linear/dunevectors.hh>

#include <dumux/assembly/coloring.hh>
#include <dumux/assembly/jacobianpattern.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/assembler.hh>

#include <dumux/parallel/multithreading.hh>
#include <dumux/parallel/parallel_for.hh>

#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Assembly
 * \brief A multi-stage (Runge-Kutta) linear system assembler (residual and Jacobian) for
 *        general discretization schemes (box, cvfe hybrid pq2/pq3, ...), dispatching to the
 *        same localDof/grid-variables local assembler as the stationary Assembler.
 * \tparam TypeTag The TypeTag
 * \tparam diffMethod The differentiation method to residual compute derivatives
 * \note Only implicit (DIRK-type) time-stepping schemes are supported: the underlying
 *       CVFELocalAssembler only has an implicit specialization.
 */
template<class TypeTag, DiffMethod diffMethod>
class MultiStageAssembler
{
    using GridGeo = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeo::GridView;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementSeed = typename GridView::Grid::template Codim<0>::EntitySeed;

    using ThisType = MultiStageAssembler<TypeTag, diffMethod>;
    using LocalAssembler = typename Detail::LocalAssemblerChooser_t<TypeTag, ThisType, diffMethod, /*isImplicit=*/true>;

public:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using StageParams = Experimental::MultiStageParams<Scalar>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using ResidualType = typename Dumux::Detail::NativeDuneVectorType<SolutionVector>::type;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridGeometry = GridGeo;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    /*!
     * \brief The constructor for instationary problems
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    MultiStageAssembler(std::shared_ptr<const Problem> problem,
                        std::shared_ptr<const GridGeometry> gridGeometry,
                        std::shared_ptr<GridVariables> gridVariables,
                        std::shared_ptr<const Experimental::MultiStageMethod<Scalar>> msMethod,
                        const SolutionVector& prevSol)
    : timeSteppingMethod_(msMethod)
    , problem_(problem)
    , gridGeometry_(gridGeometry)
    , gridVariables_(gridVariables)
    , prevSol_(&prevSol)
    {
        if (!timeSteppingMethod_->implicit())
            DUNE_THROW(Dune::NotImplemented,
                "MultiStageAssembler only supports implicit (DIRK-type) time-stepping schemes; "
                "the underlying CVFELocalAssembler has no explicit-scheme specialization.");

        enableMultithreading_ = SupportsColoring<typename GridGeometry::DiscretizationMethod>::value
            && Grid::Capabilities::supportsMultithreading(gridGeometry_->gridView())
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
        markConstrainedDofs_();

        spatialOperatorEvaluations_.back() = 0.0;
        temporalOperatorEvaluations_.back() = 0.0;

        if (stageParams_->size() != spatialOperatorEvaluations_.size())
            DUNE_THROW(Dune::InvalidStateException, "Wrong number of residuals");

        assemble_([&](const Element& element)
        {
            LocalAssembler localAssembler(*this, element, curSol);
            localAssembler.assembleJacobianAndResidual(
                *jacobian_, *residual_, *gridVariables_,
                *stageParams_,
                temporalOperatorEvaluations_.back(),
                spatialOperatorEvaluations_.back(),
                constrainedDofs_
            );
        });

        // assemble the full residual for the time integration stage
        auto constantResidualComponent = (*residual_);
        constantResidualComponent = 0.0;
        for (std::size_t k = 0; k < stageParams_->size()-1; ++k)
        {
            if (!stageParams_->skipTemporal(k))
                constantResidualComponent.axpy(stageParams_->temporalWeight(k), temporalOperatorEvaluations_[k]);
            if (!stageParams_->skipSpatial(k))
                constantResidualComponent.axpy(stageParams_->spatialWeight(k), spatialOperatorEvaluations_[k]);
        }

        // masked summation of constant residual component onto this stage's residual component
        for (std::size_t i = 0; i < constantResidualComponent.size(); ++i)
            for (std::size_t ii = 0; ii < constantResidualComponent[i].size(); ++ii)
                (*residual_)[i][ii] += constrainedDofs_[i][ii] > 0.5 ? 0.0 : constantResidualComponent[i][ii];

        applyDirichletConstraints_(curSol);
    }

    /*!
     * \brief Assembles the residual of the current stage for the given solution.
     */
    void assembleResidual(const SolutionVector& curSol)
    {
        resetResidual_();
        markConstrainedDofs_();

        if (stageParams_->size() != spatialOperatorEvaluations_.size())
            DUNE_THROW(Dune::InvalidStateException, "Wrong number of residuals");

        assembleOperatorEvaluations_(curSol);

        const auto k = stageParams_->size() - 1;
        (*residual_) = 0.0;
        residual_->axpy(stageParams_->temporalWeight(k), temporalOperatorEvaluations_.back());
        residual_->axpy(stageParams_->spatialWeight(k), spatialOperatorEvaluations_.back());

        auto constantResidualComponent = (*residual_);
        constantResidualComponent = 0.0;
        for (std::size_t i = 0; i < k; ++i)
        {
            if (!stageParams_->skipTemporal(i))
                constantResidualComponent.axpy(stageParams_->temporalWeight(i), temporalOperatorEvaluations_[i]);
            if (!stageParams_->skipSpatial(i))
                constantResidualComponent.axpy(stageParams_->spatialWeight(i), spatialOperatorEvaluations_[i]);
        }

        for (std::size_t i = 0; i < constantResidualComponent.size(); ++i)
            for (std::size_t ii = 0; ii < constantResidualComponent[i].size(); ++ii)
                (*residual_)[i][ii] += constrainedDofs_[i][ii] > 0.5 ? 0.0 : constantResidualComponent[i][ii];

        applyDirichletResidual_(curSol);
    }

    /*!
     * \brief The version without arguments uses the default constructor to create
     *        the jacobian and residual objects in this assembler if you don't need them outside this class
     */
    void setLinearSystem()
    {
        jacobian_ = std::make_shared<JacobianMatrix>();
        jacobian_->setBuildMode(JacobianMatrix::random);
        residual_ = std::make_shared<ResidualType>();

        setResidualSize_(*residual_);
        setJacobianPattern_();
    }

    /*!
     * \brief Resizes jacobian and residual and recomputes colors
     */
    void updateAfterGridAdaption()
    {
        setResidualSize_(*residual_);
        setJacobianPattern_();
        maybeComputeColors_();
    }

    //! Returns the number of degrees of freedom
    std::size_t numDofs() const
    { return gridGeometry_->numDofs(); }

    //! The problem
    const Problem& problem() const
    { return *problem_; }

    //! The global finite volume geometry
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

    //! The grid discretization
    const GridGeometry& gridDiscretization() const
    { return *gridGeometry_; }

    //! The gridview
    const GridView& gridView() const
    { return gridGeometry().gridView(); }

    //! The global grid variables
    GridVariables& gridVariables()
    { return *gridVariables_; }

    //! The global grid variables
    const GridVariables& gridVariables() const
    { return *gridVariables_; }

    //! The jacobian matrix
    JacobianMatrix& jacobian()
    { return *jacobian_; }

    //! The residual vector (rhs)
    ResidualType& residual()
    { return *residual_; }

    //! The solution of the previous time step
    const SolutionVector& prevSol() const
    { return *prevSol_; }

    /*!
     * \brief Create a local residual object (used by the local assembler)
     * \note unweighted (plain) construction; the multi-stage local assembler evaluates raw
     *       storage/flux terms and the Butcher-tableau weighting happens in this class.
     */
    LocalResidual localResidual() const
    { return LocalResidual(problem_.get(), nullptr); }

    /*!
     * \brief Update the grid variables
     */
    void updateGridVariables(const SolutionVector &cursol)
    { gridVariables().update(cursol); }

    /*!
     * \brief Reset the gridVariables
     */
    void resetTimeStep(const SolutionVector &cursol)
    {
        gridVariables().resetTimeStep(cursol);
        this->clearStages();
    }

    void clearStages()
    {
        spatialOperatorEvaluations_.clear();
        temporalOperatorEvaluations_.clear();
        stageParams_.reset();
    }

    template<class StageParamsArg>
    void prepareStage(SolutionVector& x, StageParamsArg params)
    {
        stageParams_ = std::move(params);
        const auto curStage = stageParams_->size() - 1;

        // in the first stage, also assemble the residual
        // at the previous time level (stage 0 residual)
        if (curStage == 1)
        {
            setProblemTime_(*problem_, stageParams_->timeAtStage(0));

            resetResidual_();

            assert(spatialOperatorEvaluations_.size() >= 0);
            if (spatialOperatorEvaluations_.size() == 0)
            {
                spatialOperatorEvaluations_.push_back(*residual_);
                temporalOperatorEvaluations_.push_back(*residual_);
                assembleOperatorEvaluations_(*prevSol_);
            }

            // we don't delete the first stage so it can be reused in a restarted
            // time integration step. The evaluations are only deleted
            // when explicitly requested by calling clearStages().
            // So if here the vector is non-empty, we don't need to evaluate again
            // (this should only occur if we are restarting time integration, e.g.
            // with a different time step size)
            else if (spatialOperatorEvaluations_.size() > 0)
            {
                updateGridVariables(x);
                spatialOperatorEvaluations_.resize(1);
                temporalOperatorEvaluations_.resize(1);
            }
        }

        if (spatialOperatorEvaluations_.size() != curStage)
            DUNE_THROW(Dune::InvalidStateException,
                "Invalid state. Maybe you forgot to call clearStages()");

        // The evaluations recorded while solving the previous stage belong to the last iterate
        // the solver assembled, which is not the stage solution for a solver that stops after an
        // update (e.g. Newton on the shift criterion) or assembles only once (a linear solver).
        if (curStage > 1)
        {
            setProblemTime_(*problem_, stageParams_->timeAtStage(curStage-1));
            assembleOperatorEvaluations_(x);
        }

        setProblemTime_(*problem_, stageParams_->timeAtStage(curStage));

        resetResidual_();

        // allocate memory for this stage
        spatialOperatorEvaluations_.push_back(*residual_);
        temporalOperatorEvaluations_.push_back(*residual_);
    }

    //! TODO get rid of this (called by Newton but shouldn't be necessary)
    bool isStationaryProblem() const
    { return false; }

    bool isImplicit() const
    { return timeSteppingMethod_->implicit(); }

    //! The temporal and spatial weight of the current stage
    std::pair<Scalar, Scalar> currentStageWeights() const
    {
        if (!stageParams_)
            DUNE_THROW(Dune::InvalidStateException, "No stage params set. Call prepareStage first.");
        const auto k = stageParams_->size() - 1;
        return {stageParams_->temporalWeight(k), stageParams_->spatialWeight(k)};
    }

private:
    //! Assemble the unweighted temporal and spatial operators at sol into the last stored evaluations
    void assembleOperatorEvaluations_(const SolutionVector& sol)
    {
        spatialOperatorEvaluations_.back() = 0.0;
        temporalOperatorEvaluations_.back() = 0.0;
        assemble_([&](const Element& element)
        {
            LocalAssembler localAssembler(*this, element, sol);
            localAssembler.assembleCurrentResidual(temporalOperatorEvaluations_.back(),
                                                   spatialOperatorEvaluations_.back());
        });
    }

    /*!
     * \brief Resizes the jacobian and sets the jacobian's sparsity pattern.
     */
    void setJacobianPattern_()
    {
        const auto numDofs = this->numDofs();
        jacobian_->setSize(numDofs, numDofs);
        getJacobianPattern<true>(gridGeometry()).exportIdx(*jacobian_);
    }

    //! Resizes the residual
    void setResidualSize_(ResidualType& res)
    { res.resize(numDofs()); }

    //! Computes the colors
    void maybeComputeColors_()
    {
        if (enableMultithreading_)
            elementSets_ = computeColoring(gridGeometry()).sets;
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

    // reset the Jacobian matrix to 0.0
    void resetJacobian_()
    {
        if(!jacobian_)
        {
            jacobian_ = std::make_shared<JacobianMatrix>();
            jacobian_->setBuildMode(JacobianMatrix::random);
            setJacobianPattern_();
        }

        *jacobian_ = 0.0;
    }

    //! Mark the dofs constrained via problem.constraints() so the multi-stage residual
    //! combination step (assembleJacobianAndResidual) excludes them from the accumulated
    //! historical-stage residual before the Dirichlet rows are stamped in below.
    void markConstrainedDofs_()
    {
        if constexpr (Detail::hasGlobalConstraints<Problem>())
        {
            for (const auto& constraintData : problem_->constraints())
            {
                const auto& info = constraintData.constraintInfo();
                const auto dofIdx = constraintData.dofIndex();
                for (int eqIdx = 0; eqIdx < info.size(); ++eqIdx)
                    if (info.isConstraintEquation(eqIdx))
                        constrainedDofs_[dofIdx][eqIdx] = 1.0;
            }
        }
    }

    //! Stamp Dirichlet constraint rows (problem.constraints()) into the Jacobian and residual.
    void applyDirichletConstraints_(const SolutionVector& curSol)
    {
        if constexpr (Detail::hasGlobalConstraints<Problem>())
        {
            for (const auto& constraintData : problem_->constraints())
            {
                const auto& info = constraintData.constraintInfo();
                const auto& values = constraintData.values();
                const auto dofIdx = constraintData.dofIndex();
                for (int eqIdx = 0; eqIdx < info.size(); ++eqIdx)
                {
                    if (info.isConstraintEquation(eqIdx))
                    {
                        const auto pvIdx = info.eqToPriVarIndex(eqIdx);
                        (*residual_)[dofIdx][eqIdx] = curSol[dofIdx][pvIdx] - values[pvIdx];

                        auto& row = (*jacobian_)[dofIdx];
                        for (auto col = row.begin(); col != row.end(); ++col)
                            row[col.index()][eqIdx] = 0.0;

                        (*jacobian_)[dofIdx][dofIdx][eqIdx][pvIdx] = 1.0;
                    }
                }
            }
        }
    }

    //! Set the residual of the Dirichlet constraint rows (problem.constraints()).
    void applyDirichletResidual_(const SolutionVector& curSol)
    {
        if constexpr (Detail::hasGlobalConstraints<Problem>())
        {
            for (const auto& constraintData : problem_->constraints())
            {
                const auto& info = constraintData.constraintInfo();
                const auto& values = constraintData.values();
                const auto dofIdx = constraintData.dofIndex();
                for (int eqIdx = 0; eqIdx < info.size(); ++eqIdx)
                    if (info.isConstraintEquation(eqIdx))
                    {
                        const auto pvIdx = info.eqToPriVarIndex(eqIdx);
                        (*residual_)[dofIdx][eqIdx] = curSol[dofIdx][pvIdx] - values[pvIdx];
                    }
            }
        }
    }

    /*!
     * \brief A method assembling something per element
     * \note Handles exceptions for parallel runs
     * \throws NumericalProblem on all processes if an exception is thrown during assembly
     */
    template<typename AssembleElementFunc>
    void assemble_(AssembleElementFunc&& assembleElement) const
    {
        // a state that will be checked on all processes
        bool succeeded = false;

        try
        {
            if (enableMultithreading_)
            {
                assert(elementSets_.size() > 0);

                for (const auto& elements : elementSets_)
                {
                    Dumux::parallelFor(elements.size(), [&](const std::size_t i)
                    {
                        const auto element = gridView().grid().entity(elements[i]);
                        assembleElement(element);
                    });
                }
            }
            else
                for (const auto& element : elements(gridView()))
                    assembleElement(element);

            succeeded = true;
        }
        catch (NumericalProblem &e)
        {
            std::cout << "rank " << gridView().comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n";
            succeeded = false;
        }

        if (gridView().comm().size() > 1)
            succeeded = gridView().comm().min(succeeded);

        if (!succeeded)
            DUNE_THROW(NumericalProblem, "A process did not succeed in linearizing the system");
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

    //! pointer to the problem to be solved
    std::shared_ptr<const Problem> problem_;

    //! the finite volume geometry of the grid
    std::shared_ptr<const GridGeometry> gridGeometry_;

    //! the variables container for the grid
    std::shared_ptr<GridVariables> gridVariables_;

    //! an observing pointer to the previous solution for instationary problems
    const SolutionVector* prevSol_ = nullptr;

    //! shared pointers to the jacobian matrix and residual
    std::shared_ptr<JacobianMatrix> jacobian_;
    std::shared_ptr<ResidualType> residual_;

    //! element sets for parallel assembly
    bool enableMultithreading_ = false;
    std::deque<std::vector<ElementSeed>> elementSets_;
};

} // namespace Dumux::Experimental

#endif
