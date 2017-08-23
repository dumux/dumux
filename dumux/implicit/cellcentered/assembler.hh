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
 * \brief An assembler for the global linear system
 *        for fully implicit models and cell-centered discretization schemes.
 */
#ifndef DUMUX_IMPLICIT_CC_ASSEMBLER_HH
#define DUMUX_IMPLICIT_CC_ASSEMBLER_HH

#include <dune/istl/matrixindexset.hh>

#include <dumux/common/timeloop.hh>
#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>

#include "localassembler.hh"

namespace Dumux {

namespace Properties
{
NEW_PROP_TAG(LocalAssembler);
}

template<class TypeTag, DiscretizationMethods DM>
class ImplicitAssemblerImplementation;

template<class TypeTag>
class CCImplicitAssembler;

//! TPFA and MPFA use the same assembler class
template<class TypeTag>
class ImplicitAssemblerImplementation<TypeTag, DiscretizationMethods::CCTpfa> : public CCImplicitAssembler<TypeTag>
{
    using CCImplicitAssembler<TypeTag>::CCImplicitAssembler;
};

template<class TypeTag>
class ImplicitAssemblerImplementation<TypeTag, DiscretizationMethods::CCMpfa> : public CCImplicitAssembler<TypeTag>
{
    using CCImplicitAssembler<TypeTag>::CCImplicitAssembler;
};

/*!
 * \ingroup ImplicitModel
 * \brief An assembler for the global linear system
 *        for fully implicit models and cell-centered discretization schemes.
 */
template<class TypeTag>
class CCImplicitAssembler
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    // using LocalAssembler = typename GET_PROP_TYPE(TypeTag, LocalAssembler);
    using LocalAssembler = CCImplicitLocalAssembler<TypeTag, DifferentiationMethods::numeric>;

public:
    using ResidualType = SolutionVector;

    //! The constructor for stationary problems
    CCImplicitAssembler(std::shared_ptr<const Problem> problem,
                        std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                        std::shared_ptr<GridVariables> gridVariables)
    : problem_(problem)
    , fvGridGeometry_(fvGridGeometry)
    , gridVariables_(gridVariables)
    {}

    //! The constructor for instationary problems
    CCImplicitAssembler(std::shared_ptr<const Problem> problem,
                        std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                        std::shared_ptr<GridVariables> gridVariables,
                        std::shared_ptr<TimeLoop<Scalar>> timeLoop)
    : problem_(problem)
    , fvGridGeometry_(fvGridGeometry)
    , gridVariables_(gridVariables)
    , localResidual_(timeLoop)
    {}

    /*!
     * \brief Assembles the global Jacobian of the residual
     *        and the residual for the current solution.
     */
    void assembleJacobianAndResidual(const SolutionVector& curSol)
    {
        if(!matrix_)
            DUNE_THROW(Dune::InvalidStateException,
                "Can't assemble jacobian because no jacobian matrix was set. "
                << "Use setLinearSystem() to do so!");

        if(!residual_)
            DUNE_THROW(Dune::InvalidStateException,
                "Can't assemble residual because no residual vector was set. "
                << "Use setLinearSystem() to do so!");

        resetSystem_();

        bool succeeded;
        // try assembling the global linear system
        try
        {
            // let the local assembler add the element contributions
            for (const auto element : elements(gridView()))
                LocalAssembler::assemble(*this, residual(), element, curSol);

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
        resetMatrix_();

        bool succeeded;
        // try assembling the global linear system
        try
        {
            // let the local assembler add the element contributions
            for (const auto& element : elements(gridView()))
                LocalAssembler::assemble(*this, element, curSol);

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
    void assembleResidual(const SolutionVector& curSol) const
    { assembleResidual(residual(), curSol); }

    //! compute the residuals
    void assembleResidual(ResidualType& r, const SolutionVector& curSol) const
    {
        r = 0.0;
        for (const auto& element : elements(gridView()))
        {
            if (element.partitionType() != Dune::GhostEntity)
            {
                const auto eIdx = fvGridGeometry().elementMapper().index(element);
                r[eIdx] += localResidual().eval(problem(),
                                                element,
                                                fvGridGeometry(),
                                                gridVariables(),
                                                curSol)[0];
            }
        }
    }

    //! computes the global residual
    Scalar globalResidual(const SolutionVector& curSol) const
    {
        ResidualType residual(numDofs());
        assembleResidual(residual, curSol);

        // calculate the square norm of the residual
        Scalar result2 = residual.two_norm2();
        if (gridView().comm().size() > 1)
            result2 = gridView().comm().sum(result2);

        using std::sqrt;
        return sqrt(result2);
    }

    /*!
     * \brief Tells the assembler which matrix and right hand side vector to use.
     *        This also resizes the containers to the required sizes and sets the
     *        sparsity pattern of the matrix.
     */
    void setLinearSystem(std::shared_ptr<JacobianMatrix> A,
                         std::shared_ptr<SolutionVector> r)
    {
        matrix_ = A;
        residual_ = r;

        // check and/or set the BCRS matrix's build mode
        if (matrix().buildMode() == JacobianMatrix::BuildMode::unknown)
            matrix().setBuildMode(JacobianMatrix::random);
        else if (matrix().buildMode() != JacobianMatrix::BuildMode::random)
            DUNE_THROW(Dune::NotImplemented, "Only BCRS matrices with random build mode are supported at the moment");

        allocateLinearSystem();
    }

    /*!
     * \brief The version without arguments uses the default constructor to create
     *        the jacobian and residual objects in this assembler.
     */
    void setLinearSystem()
    {
        matrix_ = std::make_shared<JacobianMatrix>();
        matrix_->setBuildMode(JacobianMatrix::random);

        residual_ = std::make_shared<SolutionVector>();

        allocateLinearSystem();
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
     * \brief Resizes the matrix and right hand side and sets the matrix' sparsity pattern.
     */
    void allocateLinearSystem()
    {
        // resize the matrix and the residual
        const auto numDofs = this->numDofs();
        matrix().setSize(numDofs, numDofs);
        residual().resize(numDofs);

        // get occupation pattern of the matrix
        Dune::MatrixIndexSet occupationPattern;
        occupationPattern.resize(numDofs, numDofs);

        for (unsigned int globalI = 0; globalI < numDofs; ++globalI)
        {
            occupationPattern.add(globalI, globalI);
            for (const auto& dataJ : fvGridGeometry().connectivityMap()[globalI])
                occupationPattern.add(dataJ.globalJ, globalI);

            // reserve index for additional user defined DOF dependencies
            // const auto& additionalDofDependencies = problem().getAdditionalDofDependencies(globalI);
            // for (auto globalJ : additionalDofDependencies)
            //     occupationPattern.add(globalI, globalJ);
        }

        // export pattern to matrix
        occupationPattern.exportIdx(matrix());
    }

    //! cell-centered schemes have one dof per cell
    std::size_t numDofs() const
    { return gridView().size(0); }

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

    JacobianMatrix& matrix()
    {
        assert(matrix_ && "The matrix has not been set yet!");
        return *matrix_;
    }

    SolutionVector& residual()
    {
        assert(residual_ && "The residual has not been set yet!");
        return *residual_;
    }

    const LocalResidual& localResidual() const
    { return localResidual_; }

private:
    // reset the global linear system of equations.
    // TODO: if partial reassemble is enabled, this means
    // that the jacobian matrix must only be erased partially!
    void resetSystem_()
    {
        residual() = 0.0;
        resetMatrix_();
    }
    void resetMatrix_()
    {
        matrix() = 0.0;
    }

    // pointer to the problem to be solved
    std::shared_ptr<const Problem> problem_;

    // the finite volume geometry of the grid
    std::shared_ptr<const FVGridGeometry> fvGridGeometry_;

    // the variables container for the grid
    std::shared_ptr<GridVariables> gridVariables_;

    // shared pointers to the jacobian matrix and residual
    std::shared_ptr<JacobianMatrix> matrix_;
    std::shared_ptr<SolutionVector> residual_;

    // class computing the residual of an element
    LocalResidual localResidual_;
};

} // namespace Dumux

#endif
