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
 * \brief Assembler class for residuals and jacobian matrices for grid-based numerical schemes.
 */
#ifndef DUMUX_ASSEMBLER_HH
#define DUMUX_ASSEMBLER_HH

#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dumux/common/timeloop.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/localview.hh>

#include "cclocalassembler.hh"
#include "boxlocalassembler.hh"
#include "felocalassembler.hh"
#include "jacobianpattern.hh"

namespace Dumux {
namespace Impl {

    template<class Assembler, DiscretizationMethod dm> struct LocalAssemblerChooser;

    // TODO: Currently the fv local residuals require type tag etc...
    // template<class Assembler> struct LocalAssemblerChooser<Assembler, DiscretizationMethod::tpfa> { using type = CCLocalAssembler<... };
    // template<class Assembler> struct LocalAssemblerChooser<Assembler, DiscretizationMethod::mpfa> { using type = CCLocalAssembler<... };
    // template<class Assembler> struct LocalAssemblerChooser<Assembler, DiscretizationMethod::box> { using type = BoxLocalAssembler<... };
    // template<class Assembler> struct LocalAssemblerChooser<Assembler, DiscretizationMethod::staggered> { using type = ... };
    template<class Assembler> struct LocalAssemblerChooser<Assembler, DiscretizationMethod::fem> { using type = FELocalAssembler<Assembler>; };

    template<class Assembler, DiscretizationMethod dm>
    using LocalAssemblerType = typename LocalAssemblerChooser<Assembler, dm>::type;

} // end namespace detail

//! Default types used for the linear system
template<class Scalar, int numEq>
struct DefaultLinearSystemTraits
{
private:
    using PrimaryVariables = Dune::FieldVector<Scalar, numEq>;
    using BlockType = Dune::FieldMatrix<Scalar, numEq, numEq>;

public:
    using Residual = Dune::BlockVector<PrimaryVariables>;
    using JacobianMatrix = Dune::BCRSMatrix<BlockType>;
};

/*!
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for grid-based numerical schemes
 * \tparam GV The grid variables type
 * \tparam LR The local residual type
 * \tparam diffMethod The differentiation method to compute derivatives
 * \tparam LST The linear system traits (types used for jacobians and residuals)
 */
template< class GV, class LR, DiffMethod diffMethod,
          class LST = DefaultLinearSystemTraits<typename GV::Scalar, GV::PrimaryVariables::size()> >
class Assembler
{
    using ThisType = Assembler<GV, LR, diffMethod, LST>;

    using GridView = typename GV::GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using TimeLoop = TimeLoopBase<typename GV::Scalar>;

public:
    //! export the types used for the linear system
    using Scalar = typename GV::Scalar;
    using JacobianMatrix = typename LST::JacobianMatrix;
    using Residual = typename LST::Residual;
    using ResidualType = Residual; // Required by NewtonSolver etc. (TODO: Which name is better!?)

    //! export the local residual type
    using LocalResidual = LR;

    //! export the problem to be solved, the geometry on which it is solved
    //! and the variables required to assemble the equations.
    using GridVariables = GV;
    using GridGeometry = typename GridVariables::GridGeometry;
    using Problem = typename GridVariables::Problem;

    //! The local assembler is the local view on this assembler class
    using LocalView = Impl::LocalAssemblerType<Assembler, GV::GridGeometry::discMethod>;

    //! Export type used for solution vectors
    // (TODO: CAN BE REMOVED IF SOLUTION VECTOR DISAPPEARS FROM INTERFACES (GRID VARS CARRY NOTION OF CURRENT SOLUTION)??)
    using SolutionVector = typename GridVariables::SolutionVector;

    /*!
     * \brief The constructor for stationary problems
     */
    Assembler(std::shared_ptr<GridVariables> gridVariables)
    : gridVariables_(gridVariables)
    , timeLoop_()
    , isStationaryProblem_(true)
    {}

    /*!
     * \brief The constructor for instationary problems
     */
    Assembler(std::shared_ptr<GridVariables> gridVariables,
              std::shared_ptr<const TimeLoop> timeLoop,
              const SolutionVector& prevSol)
    : gridVariables_(gridVariables)
    , timeLoop_(timeLoop)
    , prevSol_(&prevSol)
    , isStationaryProblem_(!timeLoop)
    {}

    /*!
     * \brief Assembles the global Jacobian of the residual and the residual around the given solution.
     * \todo TODO: The current solution is actually not used, since we assume the grid variables to carry
     *             the current state of the system. This is here for compatibility with the current implementation
     *             of the NewtonSolver!
     */
    template<class PartialReassembler = DefaultPartialReassembler>
    void assembleJacobianAndResidual(const SolutionVector& curSol,
                                     const PartialReassembler* partialReassembler = nullptr)
    {
        checkAssemblerState_();
        resetJacobian_(partialReassembler);
        resetResidual_();

        assemble_([&](const Element& element)
        {
            auto localAssembler = localView(*this);
            localAssembler.bind(element);
            localAssembler.assembleJacobianAndResidual(*jacobian_, *residual_, partialReassembler);
        });

        enforcePeriodicConstraints_(*jacobian_, *residual_, gridVariables_->gridGeometry());
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
            auto localAssembler = localView(*this);
            localAssembler.bind(element);
            localAssembler.assembleJacobianAndResidual(*jacobian_);
        });
    }

    //! compute the residuals using the internal residual
    void assembleResidual(const SolutionVector& curSol)
    {
        resetResidual_();
        assembleResidual(*residual_, curSol);
    }

    //! assemble a residual
    void assembleResidual(Residual& r, const SolutionVector& curSol) const
    {
        checkAssemblerState_();

        // update the grid variables for the case of active caching
        gridVariables_->update(curSol);

        assemble_([&](const Element& element)
        {
            auto localAssembler = localView(*this);
            localAssembler.bind(element);
            localAssembler.assembleResidual(r);
        });
    }

    //! compute the residual and return it's vector norm
    Scalar residualNorm(const SolutionVector& curSol) const
    {
        Residual residual(numDofs());
        assembleResidual(residual, curSol);

        // for box communicate the residual with the neighboring processes
        // TODO: INSTEAD OF ISBOX, DETERMINE GENERALLY IF THE SCHEME NEEDS THIS AND UNCOMMENT
        // if (isBox && gridView().comm().size() > 1)
        // {
        //     using VertexMapper = typename GridGeometry::VertexMapper;
        //     VectorCommDataHandleSum<VertexMapper, SolutionVector, GridGeometry::GridView::dimension>
        //         sumResidualHandle(gridGeometry_->vertexMapper(), residual);
        //     gridView().communicate(sumResidualHandle,
        //                            Dune::InteriorBorder_InteriorBorder_Interface,
        //                            Dune::ForwardCommunication);
        // }

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
        // TODO: HOW TO DETERMINE PATTERN DEPENDING ON TIME SCHEME??
        const auto occupationPattern = getJacobianPattern</*isImpicit*/true>(gridGeometry());

        // export pattern to jacobian
        occupationPattern.exportIdx(*jacobian_);
    }

    //! Resizes the residual
    void setResidualSize()
    { residual_->resize(numDofs()); }

    //! Returns the number of degrees of freedom
    std::size_t numDofs() const
    { return gridGeometry().numDofs(); }

    //! The problem
    const Problem& problem() const
    { return gridVariables_->problem(); }

    //! The global finite volume geometry
    const GridGeometry& gridGeometry() const
    { return gridVariables_->gridGeometry(); }

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

    //! Return the time loop
    const TimeLoop& timeLoop() const
    { assert(timeLoop_); return *timeLoop_; }

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
     * \brief Update the grid variables
     */
    void updateGridVariables(const SolutionVector& cursol)
    {
        gridVariables().update(cursol);
    }

    /*!
     * \brief Reset the gridVariables
     */
    void resetTimeStep(const SolutionVector& cursol)
    {
        gridVariables().resetTimeStep(cursol);
    }

protected:
    // reset the residual vector to 0.0
    void resetResidual_()
    {
        if (!residual_)
        {
            residual_ = std::make_shared<Residual>();
            setResidualSize();
        }

        (*residual_) = 0.0;
    }

    // reset the Jacobian matrix to 0.0
    template <class PartialReassembler = DefaultPartialReassembler>
    void resetJacobian_(const PartialReassembler *partialReassembler = nullptr)
    {
        if (!jacobian_)
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

    template<class GG>
    void enforcePeriodicConstraints_(JacobianMatrix& jac, SolutionVector& res, const GG& gridGeometry)
    { /*TODO: Implement*/ }

private:
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
