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

#include "diffmethod.hh"
#include "boxlocalassembler.hh"
#include "cclocalassembler.hh"

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for finite volume schemes
 * \tparam TypeTag the TypeTag
 * \tparam diffMethod the differentiation method to residual compute derivatives
 * \tparam isImplicit if to use an implicit or explicit time discretization
 */
template<class TypeTag, DiffMethod diffMethod, bool isImplicit = true>
class FVAssembler
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using VertexMapper = typename GET_PROP_TYPE(TypeTag, VertexMapper);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using TimeLoop = TimeLoopBase<Scalar>;

    static constexpr int dim = GridView::dimension;
    static constexpr bool isBox = GET_PROP_VALUE(TypeTag, DiscretizationMethod) == DiscretizationMethods::Box;
    using LocalAssembler = std::conditional_t<isBox, BoxLocalAssembler<TypeTag, diffMethod, isImplicit>,
                                                     CCLocalAssembler<TypeTag, diffMethod, isImplicit>>;

public:
    using ResidualType = SolutionVector;

    //! The constructor for stationary problems
    FVAssembler(std::shared_ptr<const Problem> problem,
                std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                std::shared_ptr<GridVariables> gridVariables)
    : problem_(problem)
    , fvGridGeometry_(fvGridGeometry)
    , gridVariables_(gridVariables)
    , stationary_(true)
    {
        static_assert(isImplicit, "Explicit assembler for stationary problem doesn't make sense!");
    }

    //! The constructor for instationary problems
    FVAssembler(std::shared_ptr<const Problem> problem,
                std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                std::shared_ptr<GridVariables> gridVariables,
                std::shared_ptr<TimeLoop> timeLoop)
    : problem_(problem)
    , fvGridGeometry_(fvGridGeometry)
    , gridVariables_(gridVariables)
    , localResidual_(timeLoop)
    , stationary_(false)
    {}

    /*!
     * \brief Assembles the global Jacobian of the residual
     *        and the residual for the current solution.
     */
    void assembleJacobianAndResidual(const SolutionVector& curSol)
    {
        if (!stationary_ && localResidual_.isStationary())
            DUNE_THROW(Dune::InvalidStateException, "Assembling instationary problem but previous solution was not set!");

        if(!jacobian_)
        {
            jacobian_ = std::make_shared<JacobianMatrix>();
            jacobian_->setBuildMode(JacobianMatrix::random);
            setJacobianPattern();
        }

        if(!residual_)
        {
            residual_ = std::make_shared<SolutionVector>();
            setResidualSize();
        }

        resetJacobian_();
        resetResidual_();

        bool succeeded;
        // try assembling the global linear system
        try
        {
            // let the local assembler add the element contributions
            for (const auto& element : elements(gridView()))
                LocalAssembler::assemble(*this, *jacobian_, *residual_, element, curSol);

            // if we get here, everything worked well
            succeeded = true;
            if (gridView().comm().size() > 1)
                succeeded = gridView().comm().min(succeeded);
        }
        // throw exception if a problem ocurred
        catch (NumericalProblem &e)
        {
            std::cout << "rank " << gridView().comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n";
            succeeded = false;
            if (gridView().comm().size() > 1)
                succeeded = gridView().comm().min(succeeded);
        }
        if (!succeeded)
            DUNE_THROW(NumericalProblem, "A process did not succeed in linearizing the system");
    }

    /*!
     * \brief Assembles only the global Jacobian of the residual.
     */
    void assembleJacobian(const SolutionVector& curSol)
    {
        if (!stationary_ && localResidual_.isStationary())
            DUNE_THROW(Dune::InvalidStateException, "Assembling instationary problem but previous solution was not set!");

        if(!jacobian_)
        {
            jacobian_ = std::make_shared<JacobianMatrix>();
            jacobian_->setBuildMode(JacobianMatrix::random);
            setJacobianPattern();
        }

        resetJacobian_();

        bool succeeded;
        // try assembling the global linear system
        try
        {
            // let the local assembler add the element contributions
            for (const auto& element : elements(gridView()))
                LocalAssembler::assemble(*this, *jacobian_, element, curSol);

            // if we get here, everything worked well
            succeeded = true;
            if (gridView().comm().size() > 1)
                succeeded = gridView().comm().min(succeeded);
        }
        // throw exception if a problem ocurred
        catch (NumericalProblem &e)
        {
            std::cout << "rank " << gridView().comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n";
            succeeded = false;
            if (gridView().comm().size() > 1)
                succeeded = gridView().comm().min(succeeded);
        }
        if (!succeeded)
            DUNE_THROW(NumericalProblem, "A process did not succeed in linearizing the system");
    }

    //! compute the residuals
    void assembleResidual(const SolutionVector& curSol)
    {
        if(!residual_)
        {
            residual_ = std::make_shared<SolutionVector>();
            setResidualSize();
        }

        assembleResidual(*residual_, curSol);
    }

    //! compute the residuals
    void assembleResidual(ResidualType& r, const SolutionVector& curSol) const
    {
        if (!stationary_ && localResidual_.isStationary())
            DUNE_THROW(Dune::InvalidStateException, "Assembling instationary problem but previous solution was not set!");

        // update the grid variables for the case of active caching
        gridVariables_->update(curSol);

        // let the local assembler add the element contributions
        for (const auto& element : elements(gridView()))
            LocalAssembler::assemble(*this, r, element, curSol);
    }

    //! computes the residual norm
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
     *        the jacobian and residual objects in this assembler.
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
     * \brief Sets the solution from which to start the time integration. Has to be
     *        called prior to assembly for time-dependent problems.
     */
    void setPreviousSolution(const SolutionVector& u)
    { localResidual_.setPreviousSolution(u); }

    /*!
     * \brief Return the solution that has been set as the previous one.
     */
    const SolutionVector& prevSol() const
    { return localResidual_.prevSol(); }

    /*!
     * \brief Resizes the jacobian and sets the jacobian' sparsity pattern.
     */

    void setJacobianPattern()
    {
        // resize the jacobian and the residual
        const auto numDofs = this->numDofs();
        jacobian_->setSize(numDofs, numDofs);

        // get occupation pattern of the jacobian
        Dune::MatrixIndexSet occupationPattern;
        occupationPattern.resize(numDofs, numDofs);

        // matrix pattern for implicit jacobians
        if (isImplicit)
            setImplicitJacobianPattern_(occupationPattern, numDofs);

        // matrix pattern for explicit jacobians -> diagonal matrix
        else
            for (unsigned int globalI = 0; globalI < numDofs; ++globalI)
                occupationPattern.add(globalI, globalI);

        // export pattern to jacobian
        occupationPattern.exportIdx(*jacobian_);
    }

    /*!
     * \brief Resizes the residual
     */
    void setResidualSize()
    { residual_->resize(numDofs()); }

    //! cell-centered schemes have one dof per cell
    std::size_t numDofs() const
    { return fvGridGeometry_->numDofs(); }

    const Problem& problem() const
    { return *problem_; }

    const FVGridGeometry& fvGridGeometry() const
    { return *fvGridGeometry_; }

    const GridView& gridView() const
    { return fvGridGeometry().gridView(); }

    GridVariables& gridVariables()
    { return *gridVariables_; }

    const GridVariables& gridVariables() const
    { return *gridVariables_; }

    JacobianMatrix& jacobian()
    {
       if (!residual_)
            DUNE_THROW(Dune::InvalidStateException, "No jacobian was set.");
        return *jacobian_;
    }

    SolutionVector& residual()
    {
        if (!residual_)
            DUNE_THROW(Dune::InvalidStateException, "No residual was set.");
        return *residual_;
    }

    const LocalResidual& localResidual() const
    { return localResidual_; }

private:
    // reset the residual to 0.0
    void resetResidual_()
    {
        (*residual_) = 0.0;
    }

    // reset the jacobian to 0.0
    void resetJacobian_()
    {
       (*jacobian_)  = 0.0;
    }

    //! Implicit jacobian pattern for cell-centered fv schemes
    template<typename T = TypeTag>
    std::enable_if_t<GET_PROP_VALUE(T, DiscretizationMethod) != DiscretizationMethods::Box, void>
    setImplicitJacobianPattern_(Dune::MatrixIndexSet& pattern, std::size_t numDofs)
    {
        for (unsigned int globalI = 0; globalI < numDofs; ++globalI)
        {
            pattern.add(globalI, globalI);
            for (const auto& dataJ : fvGridGeometry().connectivityMap()[globalI])
                pattern.add(dataJ.globalJ, globalI);

            // reserve index for additional user defined DOF dependencies
            // const auto& additionalDofDependencies = problem().getAdditionalDofDependencies(globalI);
            // for (auto globalJ : additionalDofDependencies)
            //     pattern.add(globalI, globalJ);
        }
    }

    //! Implicit jacobian pattern for vertex-centered fv schemes
    template<typename T = TypeTag>
    std::enable_if_t<GET_PROP_VALUE(T, DiscretizationMethod) == DiscretizationMethods::Box, void>
    setImplicitJacobianPattern_(Dune::MatrixIndexSet& pattern, std::size_t numDofs)
    {
        for (const auto& element : elements(fvGridGeometry().gridView()))
        {
            for (unsigned int vIdx = 0; vIdx < element.subEntities(dim); ++vIdx)
            {
                const auto globalI = fvGridGeometry().vertexMapper().subIndex(element, vIdx, dim);
                for (unsigned int vIdx2 = vIdx; vIdx2 < element.subEntities(dim); ++vIdx2)
                {
                    const auto globalJ = fvGridGeometry().vertexMapper().subIndex(element, vIdx2, dim);
                    pattern.add(globalI, globalJ);
                    pattern.add(globalJ, globalI);
                }
            }
        }
    }

    // pointer to the problem to be solved
    std::shared_ptr<const Problem> problem_;

    // the finite volume geometry of the grid
    std::shared_ptr<const FVGridGeometry> fvGridGeometry_;

    // the variables container for the grid
    std::shared_ptr<GridVariables> gridVariables_;

    // shared pointers to the jacobian matrix and residual
    std::shared_ptr<JacobianMatrix> jacobian_;
    std::shared_ptr<SolutionVector> residual_;

    // class computing the residual of an element
    LocalResidual localResidual_;

    // if this assembler is assembling a time dependent problem
    bool stationary_;
};

} // namespace Dumux

#endif
