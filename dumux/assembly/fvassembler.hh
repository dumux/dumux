// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief A linear system assembler (residual and Jacobian) for finite volume schemes
 */
#ifndef DUMUX_FV_ASSEMBLER_HH
#define DUMUX_FV_ASSEMBLER_HH

#include <type_traits>

#include <dune/istl/matrixindexset.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/parallel/vertexhandles.hh>

#include "jacobianpattern.hh"
#include "diffmethod.hh"
#include "boxlocalassembler.hh"
#include "cclocalassembler.hh"

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for finite volume schemes (box, tpfa, mpfa, ...)
 * \tparam TypeTag The TypeTag
 * \tparam diffMethod The differentiation method to residual compute derivatives
 * \tparam isImplicit Specifies whether the time discretization is implicit or not not (i.e. explicit)
 */
template<class TypeTag, DiffMethod diffMethod, bool isImplicit = true>
class FVAssembler
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using VertexMapper = typename GET_PROP_TYPE(TypeTag, VertexMapper);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using Element = typename GridView::template Codim<0>::Entity;
    using TimeLoop = TimeLoopBase<Scalar>;

    static constexpr bool isBox = GET_PROP_VALUE(TypeTag, DiscretizationMethod) == DiscretizationMethods::Box;

    using ThisType = FVAssembler<TypeTag, diffMethod, isImplicit>;
    using LocalAssembler = std::conditional_t<isBox, BoxLocalAssembler<TypeTag, ThisType, diffMethod, isImplicit>,
                                                     CCLocalAssembler<TypeTag, ThisType, diffMethod, isImplicit>>;

public:
    using ResidualType = SolutionVector;

    /*!
     * \brief The constructor for stationary problems
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    FVAssembler(std::shared_ptr<const Problem> problem,
                std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                std::shared_ptr<GridVariables> gridVariables)
    : problem_(problem)
    , fvGridGeometry_(fvGridGeometry)
    , gridVariables_(gridVariables)
    , timeLoop_()
    , isStationaryProblem_(true)
    {
        static_assert(isImplicit, "Explicit assembler for stationary problem doesn't make sense!");
    }

    /*!
     * \brief The constructor for instationary problems
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    FVAssembler(std::shared_ptr<const Problem> problem,
                std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                std::shared_ptr<GridVariables> gridVariables,
                std::shared_ptr<const TimeLoop> timeLoop)
    : problem_(problem)
    , fvGridGeometry_(fvGridGeometry)
    , gridVariables_(gridVariables)
    , timeLoop_(timeLoop)
    , isStationaryProblem_(false)
    {}

    /*!
     * \brief Assembles the global Jacobian of the residual
     *        and the residual for the current solution.
     */
    void assembleJacobianAndResidual(const SolutionVector& curSol)
    {
        checkAssemblerState_();
        resetJacobian_();
        resetResidual_();

        assemble_([&](const Element& element)
        {
            LocalAssembler localAssembler(*this, element, curSol);
            localAssembler.assembleJacobianAndResidual(*jacobian_, *residual_, *gridVariables_);
        });
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

        // update the grid variables for the case of active caching
        gridVariables_->update(curSol);

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

        // for box communicate the residual with the neighboring processes
        if (isBox && gridView().comm().size() > 1)
        {
            VertexHandleSum<typename SolutionVector::block_type, SolutionVector, VertexMapper>
            sumResidualHandle(residual, fvGridGeometry_->vertexMapper());
            gridView().communicate(sumResidualHandle,
                                   Dune::InteriorBorder_InteriorBorder_Interface,
                                   Dune::ForwardCommunication);
        }

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
        const auto occupationPattern = getJacobianPattern<isImplicit>(fvGridGeometry());

        // export pattern to jacobian
        occupationPattern.exportIdx(*jacobian_);
    }

    //! Resizes the residual
    void setResidualSize()
    { residual_->resize(numDofs()); }

    //! Returns the number of degrees of freedom
    std::size_t numDofs() const
    { return fvGridGeometry_->numDofs(); }

    //! The problem
    const Problem& problem() const
    { return *problem_; }

    //! The global finite volume geometry
    const FVGridGeometry& fvGridGeometry() const
    { return *fvGridGeometry_; }

    //! The gridview
    const GridView& gridView() const
    { return fvGridGeometry().gridView(); }

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

    /*!
     * \brief Set time loop for instationary problems
     * \note calling this turns this into a stationary assembler
     */
    void setTimeManager(std::shared_ptr<const TimeLoop> timeLoop)
    { timeLoop_ = timeLoop_; isStationaryProblem_ = true; }

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

    // reset the jacobian vector to 0.0
    void resetJacobian_()
    {
        if(!jacobian_)
        {
            jacobian_ = std::make_shared<JacobianMatrix>();
            jacobian_->setBuildMode(JacobianMatrix::random);
            setJacobianPattern();
        }

       (*jacobian_)  = 0.0;
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

    //! the finite volume geometry of the grid
    std::shared_ptr<const FVGridGeometry> fvGridGeometry_;

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
};

} // namespace Dumux

#endif
