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
#include <optional>
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

#include <dumux/experimental/timestepping/timelevel.hh>
#include <dumux/experimental/timestepping/timeintegrationmethod.hh>
#include <dumux/experimental/timestepping/implicitsolve.hh>

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
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using ResidualType = typename Dumux::Detail::NativeDuneVectorType<SolutionVector>::type;
    using SolveParams = Experimental::ImplicitSolveParams<Scalar, ResidualType>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridGeometry = GridGeo;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    /*!
     * \brief The constructor for instationary problems
     * \param method The time integration method (only implicit methods are supported)
     * \param prevSol The solution at the beginning of the time step (used by solvers that reset a time step)
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    MultiStageAssembler(std::shared_ptr<const Problem> problem,
                        std::shared_ptr<const GridGeometry> gridGeometry,
                        std::shared_ptr<GridVariables> gridVariables,
                        std::shared_ptr<const Experimental::TimeIntegrationMethod<Scalar>> method,
                        const SolutionVector& prevSol)
    : timeSteppingMethod_(method)
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
        if (!solveParams_)
            DUNE_THROW(Dune::InvalidStateException, "No implicit solve set. Call setImplicitSolve first.");

        resetJacobian_();
        resetResidual_();
        markConstrainedDofs_();

        assemble_([&](const Element& element)
        {
            LocalAssembler localAssembler(*this, element, curSol);
            localAssembler.assembleJacobianAndResidual(
                *jacobian_, *residual_, *gridVariables_,
                solveParams_->weights,
                constrainedDofs_
            );
        });

        addConstant_();
        applyDirichletConstraints_(curSol);
    }

    /*!
     * \brief Assembles the residual of the implicit solve for the given solution.
     */
    void assembleResidual(const SolutionVector& curSol)
    {
        if (!solveParams_)
            DUNE_THROW(Dune::InvalidStateException, "No implicit solve set. Call setImplicitSolve first.");

        resetResidual_();
        markConstrainedDofs_();

        ResidualType storage, spatial;
        assembleOperators_(curSol, {.storage = &storage, .spatial = &spatial});
        (*residual_) = 0.0;
        residual_->axpy(solveParams_->weights.storage, storage);
        residual_->axpy(solveParams_->weights.spatial, spatial);

        addConstant_();
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
     * \note unweighted (plain) construction; the local assembler evaluates the raw
     *       storage/flux terms and applies the operator weights of the implicit solve.
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
    { gridVariables().resetTimeStep(cursol); }

    /*!
     * \brief Set the weights, the time level and the contribution of known states for the next implicit solve
     */
    void setImplicitSolve(const SolutionVector&, const SolveParams& params)
    {
        if (params.weights.inertia != 0.0)
            DUNE_THROW(Dune::NotImplemented, "Inertia terms in the multi-stage assembler");

        setProblemTime_(*problem_, params.timeLevel.current());
        solveParams_ = params;
    }

    /*!
     * \brief Evaluate the unweighted operators at the given solution and time level
     * \note The grid variables have to be updated to the given solution.
     */
    void assembleOperators(const SolutionVector& sol,
                           const Experimental::TimeLevel<Scalar>& timeLevel,
                           Experimental::OperatorEvaluations<ResidualType> evaluations)
    {
        if (evaluations.inertia)
            DUNE_THROW(Dune::NotImplemented, "Inertia terms in the multi-stage assembler");

        setProblemTime_(*problem_, timeLevel.current());
        assembleOperators_(sol, evaluations);
    }

    //! TODO get rid of this (called by Newton but shouldn't be necessary)
    bool isStationaryProblem() const
    { return false; }

    bool isImplicit() const
    { return timeSteppingMethod_->implicit(); }

    //! The operator weights of the current implicit solve
    const Experimental::OperatorWeights<Scalar>& currentOperatorWeights() const
    {
        if (!solveParams_)
            DUNE_THROW(Dune::InvalidStateException, "No implicit solve set. Call setImplicitSolve first.");
        return solveParams_->weights;
    }

private:
    //! Assemble the unweighted operators at the given solution at the current problem time
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

        assemble_([&](const Element& element)
        {
            LocalAssembler localAssembler(*this, element, sol);
            localAssembler.assembleCurrentResidual(storage, spatial);
        });
    }

    //! Add the contribution of known states to all equations except the constraint equations
    void addConstant_()
    {
        if (!solveParams_->constant)
            return;

        const auto& constant = *solveParams_->constant;
        for (std::size_t i = 0; i < constant.size(); ++i)
            for (std::size_t ii = 0; ii < constant[i].size(); ++ii)
                (*residual_)[i][ii] += constrainedDofs_[i][ii] > 0.5 ? 0.0 : constant[i][ii];
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

    //! Mark the dofs constrained via problem.constraints(), which do not receive the contribution of known states
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

    std::shared_ptr<const Experimental::TimeIntegrationMethod<Scalar>> timeSteppingMethod_;
    ResidualType constrainedDofs_;
    std::optional<SolveParams> solveParams_;

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
