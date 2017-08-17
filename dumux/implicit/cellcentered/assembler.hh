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

namespace Dumux {

//! TPFA and MPFA use the same assembler class
template<class TypeTag>
class ImplicitAssemblerImplementation<TypeTag, DiscretizationMethods::CCTpfa> : CCImplicitAssembler<TypeTag>
{
    using CCImplicitAssembler<TypeTag>::CCImplicitAssembler;
};

template<class TypeTag>
class ImplicitAssemblerImplementation<TypeTag, DiscretizationMethods::CCMpfa> : CCImplicitAssembler<TypeTag>
{
    using CCImplicitAssembler<TypeTag>::CCImplicitAssembler;
};

/*!
 * \ingroup ImplicitModel
 * \brief An assembler for the global linear system
 *        for fully implicit models and cell-centered discretization schemes.
 */
template<class TypeTag>
class CCImplicitAssembler<TypeTag, DiscretizationMethods::CCTpfa>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using AssemblyMap = typename GET_PROP_TYPE(TypeTag, AssemblyMap);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using GridFvGeometry = typename GET_PROP_TYPE(TypeTag, GridFvGeometry);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

    //! the classes containing the actual assemble routines
    //! they need to be friend as well
    using NewtonAssembler = CCImplicitNewtonAssembler<TypeTag>;
    friend NewtonAssembler;

public:
    //! The constructor
    CCImplicitAssembler(std::shared_ptr<const Problem> problem,
                        std::shared_ptr<const GridFvGeometry> gridFvGeometry,
                        std::shared_ptr<GridVariables> gridVariables,
                        std::shared_ptr<JacobianMatrix> matrix,
                        std::shared_ptr<SolutionVector> rhs,
                        const bool useAssemblyMap = true)
    : problem_(problem)
    , gridFvGeometry_(gridFvGeometry)
    , gridVariables_(gridVariables)
    , matrix_(matrix)
    , rhs_(rhs)
    {
        // set the BCRS matrix's build mode
        A.setBuildMode(JacobianMatrix::random);

        // export sparsity pattern to matrix
        createMatrix();

        // maybe build assembly map
        if (useAssemblyMap)
            assemblyMap_.init(globalFvGeometry);
    }

    /*!
     * \brief Assemble the global Jacobian of the residual
     *        and the residual for the current solution.
     *
     * This method called by Newton's method. In the Newton's method,
     * the right hand side is the residual at the last time step. Also,
     * the incorporation of the BCs is different than for linear problems.
     */
    void assembleNewton()
    {
        // instantiate the class containing the assembly algorithm
        NewtonAssembler newtonAssembler;

        bool succeeded;
        // try assembling the global linear system
        try
        {
            // reset all entries to zero
            resetSystem_();

            // let the actual assembler set up the global linear system
            newtonAssembler.assemble(*this);

            // if we get here, everything worked well
            succeeded = true;
            if (gridView_().comm().size() > 1)
                succeeded = gridView_().comm().min(succeeded);
        }
        // throw exception if a problem ocurred
        catch (NumericalProblem &e)
        {
            std::cout << "rank " << problem_().gridView().comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n";
            succeeded = false;
            if (gridView_().comm().size() > 1)
                succeeded = gridView_().comm().min(succeeded);
        }
        if (!succeeded)
            DUNE_THROW(NumericalProblem, "A process did not succeed in linearizing the system");
    }

    /*!
     * \brief Set sizes and sparsity pattern of the matrix.
     *
     * This method has to be called during initialization and
     * whenever the grid has been changed.
     */
    void createMatrix()
    {
        const auto numDofs = numDofs();

        // resize the matrix and the rhs
        matrix_().setSize(numDofs, numDofs);
        rhs_().resize(numDofs);

        // get occupation pattern of the matrix
        Dune::MatrixIndexSet occupationPattern;
        occupationPattern.resize(numDofs, numDofs);

        for (unsigned int globalI = 0; globalI < numDofs; ++globalI)
        {
            occupationPattern.add(globalI, globalI);
            for (const auto& dataJ : assemblyMap_()[globalI])
                occupationPattern.add(dataJ.globalJ, globalI);

            // reserve index for additional user defined DOF dependencies
            const auto& additionalDofDependencies = problem_().getAdditionalDofDependencies(globalI);
            for (auto globalJ : additionalDofDependencies)
                occupationPattern.add(globalI, globalJ);
        }

        // export pattern to matrix
        occupationPattern.exportIdx(matrix_());

        // initialize the global system with zeros
        // TODO: Do we need to do this here?
        //       At the beginning of the assembly it is done anyway!
        resetSystem_();
    }

    //! cell-centered schemes have one dof per cell
    std::size_t numDofs() const
    { return gridView_().size(0); }

private:
    // reset the global linear system of equations.
    // TODO: if partial reassemble is enabled, this means
    // that the jacobian matrix must only be erased partially!
    void resetSystem_()
    {
        rhs_() = 0.0;
        matrix_() = 0;
    }

    /*!
     * \brief Computes the epsilon used for numeric differentiation
     *        for a given value of a primary variable.
     *
     * \param priVar The value of the primary variable
     */
    Scalar numericEpsilon(const Scalar priVar) const
    {
        // define the base epsilon as the geometric mean of 1 and the
        // resolution of the scalar type. E.g. for standard 64 bit
        // floating point values, the resolution is about 10^-16 and
        // the base epsilon is thus approximately 10^-8.
        /*
        static const Scalar baseEps
            = Dumux::geometricMean<Scalar>(std::numeric_limits<Scalar>::epsilon(), 1.0);
        */
        static const Scalar baseEps = 1e-10;
        assert(std::numeric_limits<Scalar>::epsilon()*1e4 < baseEps);
        // the epsilon value used for the numeric differentiation is
        // now scaled by the absolute value of the primary variable...
        return baseEps*(std::abs(priVar) + 1.0);
    }

    const Problem& problem_() const
    { return *problem; }

    const GridFvGeometry& gridFvGeometry_() const
    { return *gridFvGeometry_; }

    const GridView& gridView_() const
    { return gridFvGeometry_().gridView(); }

    GridVariables& gridVariables_()
    { return *gridVariables_; }

    JacobianMatrix& matrix_()
    { return *matrix_; }

    SolutionVector& rhs_()
    { return rhs_; }

    AssemblyMap& assemblyMap_()
    { return assemblyMap_; }

    LocalResidual& localResidual_()
    { return localResidual_; }

    // pointer to the problem to be solved
    std::shared_ptr<const Problem> problem_;

    // the finite volume geometry of the grid
    std::shared_ptr<const GridFvGeometry> gridFvGeometry_;

    // the variables container for the grid
    std::shared_ptr<GridVariables> gridVariables_;

    // shared pointers to the jacobian matrix and right hand side
    std::shared_ptr<JacobianMatrix> matrix_;
    std::shared_ptr<SolutionVector> rhs_;

    // assembly map needed for derivative calculations
    AssemblyMap assemblyMap_;

    // class computing the residual of an element
    LocalResidual localResidual_;
};

} // namespace Dumux

#endif
