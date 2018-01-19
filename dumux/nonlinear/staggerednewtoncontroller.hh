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
 * \ingroup Nonlinear
 * \ingroup StaggeredDiscretization
 * \brief A newton controller for staggered finite volume schemes
 */
#ifndef DUMUX_STAGGERED_NEWTON_CONTROLLER_HH
#define DUMUX_STAGGERED_NEWTON_CONTROLLER_HH

#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/nonlinear/newtoncontroller.hh>
#include <dumux/linear/linearsolveracceptsmultitypematrix.hh>
#include <dumux/linear/matrixconverter.hh>

namespace Dumux {

/*!
 * \ingroup Nonlinear
 * \ingroup StaggeredDiscretization
 * \brief A newton controller for staggered finite volume schemes
 */

template <class Scalar,
          class Comm = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> >
class StaggeredNewtonController : public NewtonController<Scalar, Comm>
{
    using ParentType = NewtonController<Scalar, Comm>;

public:
    using ParentType::ParentType;

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * If the linear solver doesn't accept multitype matrices we copy the matrix
     * into a 1x1 block BCRS matrix for solving.
     *
     * \param ls the linear solver
     * \param A The matrix of the linear system of equations
     * \param x The vector which solves the linear system
     * \param b The right hand side of the linear system
     */
    template<class LinearSolver, class JacobianMatrix, class SolutionVector>
    void solveLinearSystem(LinearSolver& ls,
                           JacobianMatrix& A,
                           SolutionVector& x,
                           SolutionVector& b)
    {
        // check matrix sizes
        assert(checkMatrix_(A) && "Sub blocks of MultiType matrix have wrong sizes!");

        try
        {
            if (this->numSteps_ == 0)
            {
                Scalar norm2 = b.two_norm2();
                if (this->comm().size() > 1)
                    norm2 = this->comm().sum(norm2);

                using std::sqrt;
                this->initialResidual_ = sqrt(norm2);
            }

            // solve by calling the appropriate implementation depending on whether the linear solver
            // is capable of handling MultiType matrices or not
            const bool converged = solveLinearSystem_(ls, A, x, b,
                                                      std::integral_constant<bool, LinearSolverAcceptsMultiTypeMatrix<LinearSolver>::value>());

            // make sure all processes converged
            int convergedRemote = converged;
            if (this->comm().size() > 1)
                convergedRemote = this->comm().min(converged);

            if (!converged) {
                DUNE_THROW(NumericalProblem,
                           "Linear solver did not converge");
            }
            else if (!convergedRemote) {
                DUNE_THROW(NumericalProblem,
                           "Linear solver did not converge on a remote process");
            }
        }
        catch (const Dune::Exception &e) {
            // make sure all processes converged
            int converged = 0;
            if (this->comm().size() > 1)
                converged = this->comm().min(converged);

            NumericalProblem p;
            p.message(e.what());
            throw p;
        }
    }



    /*!
    * \brief Update the current solution with a delta vector.
    *
    * The error estimates required for the newtonConverged() and
    * newtonProceed() methods should be updated inside this method.
    *
    * Different update strategies, such as line search and chopped
    * updates can be implemented. The default behavior is just to
    * subtract deltaU from uLastIter, i.e.
    * \f[ u^{k+1} = u^k - \Delta u^k \f]
    *
    * \param assembler The assembler for Jacobian and residual
    * \param uCurrentIter The solution vector after the current iteration
    * \param uLastIter The solution vector after the last iteration
    * \param deltaU The delta as calculated from solving the linear
    *               system of equations. This parameter also stores
    *               the updated solution.
    */
    template<class JacobianAssembler, class SolutionVector>
    void newtonUpdate(JacobianAssembler& assembler,
                      SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const SolutionVector &deltaU)
    {
        if (this->enableShiftCriterion_)
            this->newtonUpdateShift(uLastIter, deltaU);

        if (this->useLineSearch_)
        {
            this->lineSearchUpdate_(assembler, uCurrentIter, uLastIter, deltaU);
        }
        else {
            using namespace Dune::Hybrid;
            forEach(integralRange(Dune::Hybrid::size(uLastIter)), [&](const auto dofTypeIdx)
            {
                for (unsigned int i = 0; i < uLastIter[dofTypeIdx].size(); ++i)
                {
                    uCurrentIter[dofTypeIdx][i] = uLastIter[dofTypeIdx][i];
                    uCurrentIter[dofTypeIdx][i] -= deltaU[dofTypeIdx][i];
                }
            });

            if (this->enableResidualCriterion_)
            {
                this->residualNorm_ = assembler.residualNorm(uCurrentIter);
                this->reduction_ = this->residualNorm_;
                this->reduction_ /= this->initialResidual_;
            }
        }

        // update the variables class to the new solution
        assembler.gridVariables().update(uCurrentIter);
    }

     /*!
     * \brief Update the maximum relative shift of the solution compared to
     *        the previous iteration.
     *
     * \param uLastIter The current iterative solution
     * \param deltaU The difference between the current and the next solution
     */
    template<class SolutionVector>
    void newtonUpdateShift(const SolutionVector &uLastIter,
                           const SolutionVector &deltaU)
    {
        this->shift_ = 0;

        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(uLastIter)), [&](const auto dofTypeIdx)
        {
            for (int i = 0; i < int(uLastIter[dofTypeIdx].size()); ++i)
            {
                auto uNewI = uLastIter[dofTypeIdx][i];
                uNewI -= deltaU[dofTypeIdx][i];

                Scalar shiftAtDof = this->relativeShiftAtDof_(uLastIter[dofTypeIdx][i], uNewI);
                this->shift_ = std::max(this->shift_, shiftAtDof);
            }
        });

        if (this->comm().size() > 1)
            this->shift_ = this->comm().max(this->shift_);
    }

private:
    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * Specialization for linear solvers that can handle MultiType matrices.
     *
     */
    template<class LinearSolver, class JacobianMatrix, class SolutionVector>
    bool solveLinearSystem_(LinearSolver& ls,
                            JacobianMatrix& A,
                            SolutionVector& x,
                            SolutionVector& b,
                            std::true_type)
    {
        // TODO: automatically derive the precondBlockLevel
        return ls.template solve</*precondBlockLevel=*/2>(A, x, b);
    }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * Specialization for linear solvers that cannot handle MultiType matrices.
     * We copy the matrix into a 1x1 block BCRS matrix before solving.
     *
     */
    template<class LinearSolver, class JacobianMatrix, class SolutionVector>
    bool solveLinearSystem_(LinearSolver& ls,
                            JacobianMatrix& A,
                            SolutionVector& x,
                            SolutionVector& b,
                            std::false_type)
    {
        // create the bcrs matrix the IterativeSolver backend can handle
        const auto M = MatrixConverter<JacobianMatrix>::multiTypeToBCRSMatrix(A);

        // get the new matrix sizes
        const std::size_t numRows = M.N();
        assert(numRows == M.M());

        // create the vector the IterativeSolver backend can handle
        const auto bTmp = VectorConverter<SolutionVector>::multiTypeToBlockVector(b);
        assert(bTmp.size() == numRows);

        // create a blockvector to which the linear solver writes the solution
        using VectorBlock = typename Dune::FieldVector<Scalar, 1>;
        using BlockVector = typename Dune::BlockVector<VectorBlock>;
        BlockVector y(numRows);

        // solve
        const bool converged = ls.solve(M, y, bTmp);

        // copy back the result y into x
        if(converged)
            VectorConverter<SolutionVector>::retrieveValues(x, y);

        return converged;
    }

    //! helper method to assure the MultiType matrix's sub blocks have the correct sizes
    template<class JacobianMatrix>
    bool checkMatrix_(const JacobianMatrix& A)
    {
        bool matrixHasCorrectSize = true;
        using namespace Dune::Hybrid;
        using namespace Dune::Indices;
        forEach(A, [&matrixHasCorrectSize](const auto& rowOfMultiTypeMatrix)
        {
            const auto numRowsLeftMostBlock = rowOfMultiTypeMatrix[_0].N();

            forEach(rowOfMultiTypeMatrix, [&matrixHasCorrectSize, &numRowsLeftMostBlock](const auto& subBlock)
            {
                if (subBlock.N() != numRowsLeftMostBlock)
                    matrixHasCorrectSize = false;
            });
        });
        return matrixHasCorrectSize;
    }

};

} // end namespace Dumux

#endif
