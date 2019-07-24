// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for finite element schemes
 */
#ifndef DUMUX_FE_ASSEMBLER_HH
#define DUMUX_FE_ASSEMBLER_HH

#include <type_traits>

#include <dune/istl/matrixindexset.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/discretization/method.hh>
#include <dumux/parallel/vertexhandles.hh>

#include "jacobianpattern.hh"
#include "diffmethod.hh"
#include "felocalassembler.hh"

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for finite element schemes
 * \tparam TypeTag The TypeTag
 * \tparamm TSB The basis of the trial space
 * \tparam diffMethod The differentiation method to residual compute derivatives
 * \tparam isImplicit Specifies whether the time discretization is implicit or not not (i.e. explicit)
 */
template<class TypeTag,
         class TSB = typename GetPropType<TypeTag, Properties::GridGeometry>::FEBasis,
         DiffMethod diffMethod = DiffMethod::numeric,
         bool isImplicit = true>
class FEAssembler
{
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using AnsatzSpaceBasis = typename GG::FEBasis;
    using GridView = typename GG::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using TimeLoop = TimeLoopBase<GetPropType<TypeTag, Properties::Scalar>>;

    using ThisType = FEAssembler<TypeTag, TSB, diffMethod, isImplicit>;
    using LocalAssembler = FELocalAssembler<TypeTag, ThisType, diffMethod, isImplicit>;

public:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridGeometry = GG;
    using TrialSpaceBasis = TSB;

    using ResidualType = SolutionVector;

    /*!
     * \brief The constructor for stationary problems
     */
    FEAssembler(std::shared_ptr<const Problem> problem,
                std::shared_ptr<const GridGeometry> gridGeometry,
                std::shared_ptr<GridVariables> gridVariables)
    : problem_(problem)
    , gridGeometry_(gridGeometry)
    , gridVariables_(gridVariables)
    , timeLoop_()
    , isStationaryProblem_(true)
    {
        static_assert(isImplicit, "Explicit assembler for stationary problem doesn't make sense!");
    }

    /*!
     * \brief The constructor for instationary problems
     */
    FEAssembler(std::shared_ptr<const Problem> problem,
                std::shared_ptr<const GridGeometry> gridGeometry,
                std::shared_ptr<GridVariables> gridVariables,
                std::shared_ptr<const TimeLoop> timeLoop)
    : problem_(problem)
    , gridGeometry_(gridGeometry)
    , gridVariables_(gridVariables)
    , timeLoop_(timeLoop)
    , isStationaryProblem_(!timeLoop)
    {}

    /*!
     * \brief Returns true if trial and ansatz space are identical.
     */
    static constexpr bool isStandardGalerkin()
    { return std::is_same<AnsatzSpaceBasis, TrialSpaceBasis>::value; }

    /*!
     * \brief Assembles the global Jacobian of the residual
     *        and the residual for the current solution.
     */
    template<class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianAndResidual(const SolutionVector& curSol, const PartialReassembler* partialReassembler = nullptr)
    {
        checkAssemblerState_();
        resetJacobian_(partialReassembler);
        resetResidual_();

        assemble_([&](const Element& element)
        {
            LocalAssembler localAssembler(*this, element, curSol);
            localAssembler.assembleJacobianAndResidual(*jacobian_, *residual_, *gridVariables_, partialReassembler);
        });

        // TODO: periodic boundaries
    }

    /*!
     * \brief Assembles only the global Jacobian of the residual.
     */
    void assembleJacobian(const SolutionVector& curSol)
    {
        checkAssemblerState_();
        resetJacobian_();

        assemble_([&](const Element& element)
        {
            LocalAssembler localAssembler(*this, element, curSol);
            localAssembler.assembleJacobianAndResidual(*jacobian_, *gridVariables_);
        });
    }

    //! compute the residuals using the internal residual
    void assembleResidual(const SolutionVector& curSol)
    {
        resetResidual_();
        assembleResidual(*residual_, curSol);
    }

    //! assemble a residual r
    void assembleResidual(ResidualType& r, const SolutionVector& curSol) const
    {
        checkAssemblerState_();

        assemble_([&](const Element& element)
        {
            LocalAssembler localAssembler(*this, element, curSol);
            localAssembler.assembleResidual(r);
        });
    }

    //! compute the residual and return it's vector norm
    Scalar residualNorm(const SolutionVector& curSol) const
    {
        ResidualType residual(numDofs());
        assembleResidual(residual, curSol);

        // TODO: implement for parallel runs
        // We'd have to to something like for box with the vertex
        // handles but also have edge handles for dofs on the edges
        if (gridView().comm().size() > 1)
            DUNE_THROW(Dune::NotImplemented, "FEM parallel residual computation");

        // calculate the square norm of the residual
        Scalar result2 = residual.two_norm2();
        if (gridView().comm().size() > 1)
            result2 = gridView().comm().sum(result2);

        using std::sqrt;
        return sqrt(result2);
    }

    /*!
     * \brief Tells the assembler which jacobian and residual to use.
     *        This also resizes the containers to the required sizes and sets the
     *        sparsity pattern of the jacobian matrix.
     */
    void setLinearSystem(std::shared_ptr<JacobianMatrix> A,
                         std::shared_ptr<SolutionVector> r)
    {
        jacobian_ = A;
        residual_ = r;

        // check and/or set the BCRS matrix's build mode
        if (jacobian_->buildMode() == JacobianMatrix::BuildMode::unknown)
            jacobian_->setBuildMode(JacobianMatrix::random);
        else if (jacobian_->buildMode() != JacobianMatrix::BuildMode::random)
            DUNE_THROW(Dune::NotImplemented, "Only BCRS matrices with random build mode are supported at the moment");

        setJacobianPattern();
        setResidualSize();
    }

    /*!
     * \brief The version without arguments uses the default constructor to create
     *        the jacobian and residual objects in this assembler if you don't need them outside this class
     */
    void setLinearSystem()
    {
        jacobian_ = std::make_shared<JacobianMatrix>();
        jacobian_->setBuildMode(JacobianMatrix::random);
        residual_ = std::make_shared<SolutionVector>();

        setJacobianPattern();
        setResidualSize();
    }

    /*!
     * \brief Resizes the jacobian and sets the jacobian' sparsity pattern.
     */
    void setJacobianPattern()
    {
        // resize the jacobian and the residual
        const auto numDofs = this->numDofs();
        jacobian_->setSize(numDofs, numDofs);

        // create occupation pattern of the jacobian
        const auto occupationPattern = getJacobianPattern<isImplicit>(gridGeometry());

        // export pattern to jacobian
        occupationPattern.exportIdx(*jacobian_);
    }

    //! Resizes the residual
    void setResidualSize()
    { residual_->resize(numDofs()); }

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
    SolutionVector& residual()
    { return *residual_; }

    //! The solution of the previous time step
    const SolutionVector& prevSol() const
    { return *prevSol_; }

    //! The basis of the trial space (standard galerkin case)
    template< bool isSG = isStandardGalerkin(),
              std::enable_if_t<isSG, int> = 0 >
    const TrialSpaceBasis& trialSpaceBasis() const
    { return gridGeometry().feBasis(); }

    //! The basis of the trial space
    template< bool isSG = isStandardGalerkin(),
              std::enable_if_t<!isSG, int> = 0 >
    const TrialSpaceBasis& trialSpaceBasis() const
    {
        if (!trialSpaceBasis_)
            DUNE_THROW(Dune::InvalidStateException, "No trial space basis set!");
        return *trialSpaceBasis_;
    }

    /*!
     * \brief Set time loop for instationary problems
     * \note calling this turns this into a stationary assembler
     */
    void setTimeLoop(std::shared_ptr<const TimeLoop> timeLoop)
    { timeLoop_ = timeLoop; isStationaryProblem_ = !static_cast<bool>(timeLoop); }

    /*!
     * \brief Sets the solution from which to start the time integration. Has to be
     *        called prior to assembly for time-dependent problems.
     */
    void setPreviousSolution(const SolutionVector& u)
    { prevSol_ = &u;  }

    /*!
     * \brief Whether we are assembling a stationary or instationary problem
     */
    bool isStationaryProblem() const
    { return isStationaryProblem_; }

    /*!
     * \brief Create a local residual object (used by the local assembler)
     */
    LocalResidual localResidual() const
    { return LocalResidual(problem_.get(), timeLoop_.get()); }

    /*!
     * \brief Update the grid variables
     */
    void updateGridVariables(const SolutionVector &cursol)
    {
        gridVariables().update(cursol);
    }

    /*!
     * \brief Reset the gridVariables
     */
    void resetTimeStep(const SolutionVector &cursol)
    {
        gridVariables().resetTimeStep(cursol);
    }

private:
    // reset the residual vector to 0.0
    void resetResidual_()
    {
        if(!residual_)
        {
            residual_ = std::make_shared<SolutionVector>();
            setResidualSize();
        }

        (*residual_) = 0.0;
    }

    // reset the Jacobian matrix to 0.0
    template <class PartialReassembler = DefaultPartialReassembler>
    void resetJacobian_(const PartialReassembler *partialReassembler = nullptr)
    {
        if(!jacobian_)
        {
            jacobian_ = std::make_shared<JacobianMatrix>();
            jacobian_->setBuildMode(JacobianMatrix::random);
            setJacobianPattern();
        }

        if (partialReassembler)
            partialReassembler->resetJacobian(*this);
        else
            *jacobian_ = 0.0;
    }

    // check if the assembler is in a correct state for assembly
    void checkAssemblerState_() const
    {
        if (!isStationaryProblem_ && !prevSol_)
            DUNE_THROW(Dune::InvalidStateException, "Assembling instationary problem but previous solution was not set!");
    }

    /*!
     * \brief A method assembling something per element
     * \note Handles exceptions for parallel runs
     * \throws NumericalProblem on all processes if something throwed during assembly
     */
    template<typename AssembleElementFunc>
    void assemble_(AssembleElementFunc&& assembleElement) const
    {
        // a state that will be checked on all processes
        bool succeeded = false;

        // try assembling using the local assembly function
        try
        {
            // let the local assembler add the element contributions
            for (const auto& element : elements(gridView()))
                assembleElement(element);

            // if we get here, everything worked well on this process
            succeeded = true;
        }
        // throw exception if a problem ocurred
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

    //! pointer to the problem to be solved
    std::shared_ptr<const Problem> problem_;

    //! the geometry of the grid
    std::shared_ptr<const GridGeometry> gridGeometry_;

    //! the variables container for the grid
    std::shared_ptr<GridVariables> gridVariables_;

    //! the time loop for instationary problem assembly
    std::shared_ptr<const TimeLoop> timeLoop_;

    //! an observing pointer to the previous solution for instationary problems
    const SolutionVector* prevSol_ = nullptr;

    //! if this assembler is assembling an instationary problem
    bool isStationaryProblem_;

    //! shared pointers to the jacobian matrix and residual
    std::shared_ptr<JacobianMatrix> jacobian_;
    std::shared_ptr<SolutionVector> residual_;

    //! shared pointer to the trial space basis
    std::shared_ptr<TrialSpaceBasis> trialSpaceBasis_;
};

} // namespace Dumux

#endif
