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

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>

#include "localassembler.hh"

namespace Dumux {

namespace Properties
{
NEW_PROP_TAG(AssemblyMap);
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
    using AssemblyMap = typename GET_PROP_TYPE(TypeTag, AssemblyMap);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    // using LocalAssembler = typename GET_PROP_TYPE(TypeTag, LocalAssembler);
     using LocalAssembler = CCImplicitLocalAssembler<TypeTag, DifferentiationMethods::numeric>;

public:
    using ResidualType = SolutionVector;

    //! The constructor
    CCImplicitAssembler(std::shared_ptr<const Problem> problem,
                        std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                        std::shared_ptr<const SolutionVector> prevSol,
                        std::shared_ptr<const SolutionVector> curSol,
                        std::shared_ptr<GridVariables> gridVariables)
    : problem_(problem)
    , fvGridGeometry_(fvGridGeometry)
    , prevSol_(prevSol)
    , curSol_(curSol)
    , gridVariables_(gridVariables)
    {
        // build assembly map
        assemblyMap_.init(*fvGridGeometry);
    }

    /*!
     * \brief Assembles the global Jacobian of the residual
     *        and the residual for the current solution.
     */
    void assembleJacobianAndResidual()
    {
        resetSystem_();

        bool succeeded;
        // try assembling the global linear system
        try
        {
            // let the local assembler add the element contributions
            for (const auto element : elements(gridView()))
                LocalAssembler::assemble(*this, residual(), element);

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
    void assembleJacobian()
    {
        resetMatrix_();

        bool succeeded;
        // try assembling the global linear system
        try
        {
            // let the local assembler add the element contributions
            for (const auto& element : elements(gridView()))
                LocalAssembler::assemble(*this, element);

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
     * \brief Tells the assembler which matrix and right hand side vector to use.
     *        This also resizes the containers to the required sizes and sets the
     *        sparsity pattern of the matrix.
     */
    void setLinearSystem(std::shared_ptr<JacobianMatrix>& A,
                         std::shared_ptr<SolutionVector>& r)
    {
        matrix_ = A;
        residual_ = r;

        // check and/or set the BCRS matrix's build mode
        if (matrix().buildMode() == JacobianMatrix::BuildMode::unknown)
            matrix().setBuildMode(JacobianMatrix::random);
        else if (matrix().buildMode() != JacobianMatrix::BuildMode::random)
            DUNE_THROW(Dune::NotImplemented, "Only BCRS matrices with random build mode are supported at the moment");

        updateLinearSystem();
    }

    /*!
     * \brief Resizes the matrix and right hand side and sets the matrix' sparsity pattern.
     */
    void updateLinearSystem()
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
            for (const auto& dataJ : assemblyMap()[globalI])
                occupationPattern.add(dataJ.globalJ, globalI);

            // reserve index for additional user defined DOF dependencies
            // const auto& additionalDofDependencies = problem().getAdditionalDofDependencies(globalI);
            // for (auto globalJ : additionalDofDependencies)
            //     occupationPattern.add(globalI, globalJ);
        }

        // export pattern to matrix
        occupationPattern.exportIdx(matrix());
    }

    //! compute the residuals
    void assembleResidual(ResidualType& r)
    {
        r = 0.0;
        for (const auto& element : elements(gridView()))
        {
            if (element.partitionType != Dune::GhostEntity)
            {
                const auto eIdx = fvGridGeometry().elementMapper().index(element);
                r[eIdx] += localResidual().eval(element,
                                                problem(),
                                                fvGridGeometry(),
                                                gridVariables(),
                                                curSol(),
                                                prevSol());
            }
        }
    }

    //! computes the global residual
    Scalar globalResidual()
    {
        ResidualType residual;
        assembleResidual(residual);

        // calculate the square norm of the residual
        Scalar result2 = residual.two_norm2();
        if (gridView().comm().size() > 1)
            result2 = gridView().comm().sum(result2);

        using std::sqrt;
        return sqrt(result2);
    }

    const Problem& problem() const
    { return *problem_; }

    const FVGridGeometry& fvGridGeometry() const
    { return *fvGridGeometry_; }

    const GridView& gridView() const
    { return fvGridGeometry().gridView(); }

    const SolutionVector& prevSol() const
    { return *prevSol_; }

    const SolutionVector& curSol() const
    { return *curSol_; }

    GridVariables& gridVariables()
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

    const AssemblyMap& assemblyMap() const
    { return assemblyMap_; }

    LocalResidual& localResidual()
    { return localResidual_; }

    //! cell-centered schemes have one dof per cell
    std::size_t numDofs() const
    { return gridView().size(0); }

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

    // previous and current solution to the problem
    std::shared_ptr<const SolutionVector> prevSol_;
    std::shared_ptr<const SolutionVector> curSol_;

    // the variables container for the grid
    std::shared_ptr<GridVariables> gridVariables_;

    // shared pointers to the jacobian matrix and residual
    std::shared_ptr<JacobianMatrix> matrix_;
    std::shared_ptr<SolutionVector> residual_;

    // assembly map needed for derivative calculations
    AssemblyMap assemblyMap_;

    // class computing the residual of an element
    LocalResidual localResidual_;
};

} // namespace Dumux

#endif
