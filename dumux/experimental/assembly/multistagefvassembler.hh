// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for finite volume schemes
 */
#ifndef DUMUX_EXPERIMENTAL_MULTISTAGE_FV_ASSEMBLER_HH
#define DUMUX_EXPERIMENTAL_MULTISTAGE_FV_ASSEMBLER_HH

#include <vector>
#include <deque>
#include <type_traits>
#include <memory>
#include <optional>

#include <dune/istl/matrixindexset.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/gridcapabilities.hh>
#include <dumux/common/typetraits/vector.hh>

#include <dumux/discretization/method.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/dunevectors.hh>

#include <dumux/assembly/coloring.hh>
#include <dumux/assembly/jacobianpattern.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/parallel/multithreading.hh>
#include <dumux/parallel/parallel_for.hh>

#include <dumux/experimental/assembly/cvfelocalassembler.hh>
#include <dumux/experimental/assembly/cclocalassembler.hh>
#include <dumux/experimental/assembly/multistagefvlocaloperator.hh>

#include <dumux/experimental/timestepping/timelevel.hh>
#include <dumux/experimental/timestepping/timeintegrationmethod.hh>
#include <dumux/experimental/timestepping/implicitsolve.hh>

namespace Dumux::Experimental::Detail {

template<class DiscretizationMethod>
struct LocalAssemblerChooser;

template<class DM>
struct LocalAssemblerChooser<DiscretizationMethods::CVFE<DM>>
{
    template<class TypeTag, class Impl, DiffMethod diffMethod>
    using type = Experimental::CVFELocalAssembler<TypeTag, Impl, diffMethod>;
};

template<>
struct LocalAssemblerChooser<DiscretizationMethods::CCMpfa>
{
    template<class TypeTag, class Impl, DiffMethod diffMethod>
    using type = Experimental::CCLocalAssembler<TypeTag, Impl, diffMethod>;
};

template<>
struct LocalAssemblerChooser<DiscretizationMethods::CCTpfa>
{
    template<class TypeTag, class Impl, DiffMethod diffMethod>
    using type = Experimental::CCLocalAssembler<TypeTag, Impl, diffMethod>;
};

template<class TypeTag, class Impl, DiffMethod diffMethod>
using LocalAssemblerChooser_t = typename LocalAssemblerChooser<
    typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod
>::template type<TypeTag, Impl, diffMethod>;

} // end namespace Dumux::Detail

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for finite volume schemes (box, tpfa, mpfa, ...)
 *        and the implicit solves of a time integration method (see ImplicitSolveParams)
 * \tparam TypeTag The TypeTag
 * \tparam diffMethod The differentiation method to residual compute derivatives
 */
template<class TypeTag, DiffMethod diffMethod>
class MultiStageFVAssembler
{
    using GridGeo = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeo::GridView;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementSeed = typename GridView::Grid::template Codim<0>::EntitySeed;

    static constexpr bool isBox = GridGeo::discMethod == DiscretizationMethods::box;

    using ThisType = MultiStageFVAssembler<TypeTag, diffMethod>;
    using LocalAssembler = typename Detail::LocalAssemblerChooser_t<TypeTag, ThisType, diffMethod>;

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
     * \param method The time integration method, which determines the Jacobian pattern
     * \param prevSol The solution at the beginning of the time step (used by solvers that reset a time step)
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    MultiStageFVAssembler(std::shared_ptr<const Problem> problem,
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
    template<class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianAndResidual(const SolutionVector& curSol, const PartialReassembler* partialReassembler = nullptr)
    {
        if (!solveParams_)
            DUNE_THROW(Dune::InvalidStateException, "No implicit solve set. Call setImplicitSolve first.");

        resetJacobian_(partialReassembler);
        resetResidual_();

        assemble_([&](const Element& element)
        {
            LocalAssembler localAssembler(*this, element, curSol);
            localAssembler.assembleJacobianAndResidual(
                *jacobian_, *residual_, *gridVariables_,
                solveParams_->weights,
                constrainedDofs_,
                partialReassembler
            );
        });

        // the contribution of known states does not enter the constraint equations
        if (solveParams_->constant)
        {
            const auto& constant = *solveParams_->constant;
            for (std::size_t i = 0; i < constant.size(); ++i)
                for (std::size_t ii = 0; ii < constant[i].size(); ++ii)
                    (*residual_)[i][ii] += constrainedDofs_[i][ii] > 0.5 ? 0.0 : constant[i][ii];
        }
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
     */
    MultiStageFVLocalOperator<LocalResidual> localResidual() const
    { return { LocalResidual{problem_.get(), nullptr} }; }

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
            DUNE_THROW(Dune::NotImplemented, "Inertia terms in the multi-stage FV assembler");

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
            DUNE_THROW(Dune::NotImplemented, "Inertia terms in the multi-stage FV assembler");

        setProblemTime_(*problem_, timeLevel.current());

        ResidualType unused;
        auto& storage = evaluations.storage ? *evaluations.storage : unused;
        auto& spatial = evaluations.spatial ? *evaluations.spatial : unused;
        for (auto* r : {&storage, &spatial})
        {
            setResidualSize_(*r);
            (*r) = 0.0;
        }

        assemble_([&](const auto& element)
        {
            LocalAssembler localAssembler(*this, element, sol);
            localAssembler.localResidual().spatialWeight(1.0, !evaluations.spatial);
            localAssembler.localResidual().temporalWeight(1.0, !evaluations.storage);
            localAssembler.assembleCurrentResidual(spatial, storage);
        });
    }

    //! TODO get rid of this (called by Newton but shouldn't be necessary)
    bool isStationaryProblem() const
    { return false; }

    bool isImplicit() const
    { return timeSteppingMethod_->implicit(); }

private:
    /*!
     * \brief Resizes the jacobian and sets the jacobian' sparsity pattern.
     */
    void setJacobianPattern_()
    {
        // resize the jacobian and the residual
        const auto numDofs = this->numDofs();
        jacobian_->setSize(numDofs, numDofs);

        // create occupation pattern of the jacobian
        if (timeSteppingMethod_->implicit())
            getJacobianPattern<true>(gridGeometry()).exportIdx(*jacobian_);
        else
            getJacobianPattern<false>(gridGeometry()).exportIdx(*jacobian_);
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
    template <class PartialReassembler = DefaultPartialReassembler>
    void resetJacobian_(const PartialReassembler *partialReassembler = nullptr)
    {
        if(!jacobian_)
        {
            jacobian_ = std::make_shared<JacobianMatrix>();
            jacobian_->setBuildMode(JacobianMatrix::random);
            setJacobianPattern_();
        }

        if (partialReassembler)
            partialReassembler->resetJacobian(*this);
        else
            *jacobian_ = 0.0;
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

        // try assembling using the local assembly function
        try
        {
            if (enableMultithreading_)
            {
                assert(elementSets_.size() > 0);

                // make this element loop run in parallel
                // for this we have to color the elements so that we don't get
                // race conditions when writing into the global matrix
                // each color can be assembled using multiple threads
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

            // if we get here, everything worked well on this process
            succeeded = true;
        }
        // throw exception if a problem occurred
        catch (NumericalProblem &e)
        {
            std::cout << "rank " << gridView().comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n";
            succeeded = false;
        }

        // make sure everything worked well on all processes
        if (gridView().comm().size() > 1)
            succeeded = gridView().comm().min(succeeded);

        // if not succeeded rethrow the error on all processes
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

} // namespace Dumux

#endif
