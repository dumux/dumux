// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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
 * \ingroup MultiDomain
 * \brief A Newton controller for mixed dimension simulations
 */
#ifndef DUMUX_MULTIDOMAIN_NEWTON_CONTROLLER_HH
#define DUMUX_MULTIDOMAIN_NEWTON_CONTROLLER_HH

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/math.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/matrixconverter.hh>
#include <dumux/nonlinear/newtoncontroller.hh>
#include <dumux/linear/linearsolveracceptsmultitypematrix.hh>

namespace Dumux {

/*!
 * \ingroup Nonlinear
 * \brief A Newton controller for mixed dimension simulations
 */
template <class Scalar, class CouplingManager,
          class Comm = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator>>
class MultiDomainNewtonController
: public NewtonController<Scalar, Comm>
{
    using ParentType = NewtonController<Scalar, Comm>;
public:
    /*!
     * \brief Constructor for stationary problems
     */
    MultiDomainNewtonController(std::shared_ptr<CouplingManager> couplingManager,
                                const Comm& comm = Dune::MPIHelper::getCollectiveCommunication(),
                                const std::string& paramGroup = "")
    : ParentType(comm, paramGroup)
    , couplingManager_(couplingManager)
    {}

    /*!
     * \brief Constructor for instationary problems
     */
    MultiDomainNewtonController(std::shared_ptr<TimeLoop<Scalar>> timeLoop,
                                std::shared_ptr<CouplingManager> couplingManager,
                                const Comm& comm = Dune::MPIHelper::getCollectiveCommunication(),
                                const std::string& paramGroup = "")
    : ParentType(timeLoop, comm, paramGroup)
    , couplingManager_(couplingManager)
    {}

    /*!
     * \brief Indicates the beginning of a Newton iteration.
     */
    template<class SolutionVector>
    void newtonBeginStep(const SolutionVector& uCurrentIter)
    {
        ParentType::newtonBeginStep(uCurrentIter);
        // update the coupling manager's solution vector
        couplingManager_->updateSolution(uCurrentIter);
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
        forEach(integralRange(Dune::Hybrid::size(uLastIter)), [&](const auto domainId)
        {
            for (int i = 0; i < int(uLastIter[domainId].size()); ++i)
            {
                auto uNewI = uLastIter[domainId][i];
                uNewI -= deltaU[domainId][i];

                Scalar shiftAtDof = this->relativeShiftAtDof_(uLastIter[domainId][i], uNewI);
                using std::max;
                this->shift_ = max(this->shift_, shiftAtDof);
            }
        });

        if (this->comm().size() > 1)
            this->shift_ = this->comm().max(this->shift_);
    }

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
     * \param assembler The assembler (needed for global residual evaluation)
     * \param uCurrentIter The solution vector after the current iteration
     * \param uLastIter The solution vector after the last iteration
     * \param deltaU The delta as calculated from solving the linear
     *               system of equations. This parameter also stores
     *               the updated solution.
     * TODO: introducte update GridVariables for normal newton controller
     * then this overload can be removed
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
            uCurrentIter = uLastIter;
            uCurrentIter -= deltaU;

            if (this->enableResidualCriterion_)
            {
                this->residualNorm_ = assembler.residualNorm(uCurrentIter);
                this->reduction_ = this->residualNorm_;
                this->reduction_ /= this->initialResidual_;
            }
            else
            {
                // If we get here, the convergence criterion does not require
                // additional residual evalutions. Thus, the grid variables have
                // not yet been updated to the new uCurrentIter.
                assembler.updateGridVariables(uCurrentIter);
            }
        }
    }

    /*!
     * \brief Indicates that one Newton iteration was finished.
     *
     * \param assembler The jacobian assembler
     * \param uCurrentIter The solution after the current Newton iteration
     * \param uLastIter The solution at the beginning of the current Newton iteration
     */
    template<class JacobianAssembler, class SolutionVector>
    void newtonEndStep(JacobianAssembler& assembler,
                       SolutionVector &uCurrentIter,
                       const SolutionVector &uLastIter)
    {
        couplingManager_->updateSolution(uCurrentIter);
        ParentType::newtonEndStep(assembler, uCurrentIter, uLastIter);
    }

    /*!
     * \brief Called if the Newton method broke down.
     * This method is called _after_ newtonEnd()
     * TODO: Fix resetTimeStep for other models
     */
    template<class Assembler, class SolutionVector>
    void newtonFail(Assembler& assembler, SolutionVector& u)
    {
        if (!assembler.isStationaryProblem())
        {
            // set solution to previous solution
            u = assembler.prevSol();

            // reset the grid variables to the previous solution
            assembler.resetTimeStep(u);

            if (this->verbose())
            {
                std::cout << "Newton solver did not converge with dt = "
                          << this->timeLoop_->timeStepSize() << " seconds. Retrying with time step of "
                          << this->timeLoop_->timeStepSize()/2.0 << " seconds\n";
            }

            // try again with dt = dt/2
            this->timeLoop_->setTimeStepSize(this->timeLoop_->timeStepSize()/2.0);
        }
        else
            DUNE_THROW(Dune::MathError, "Newton solver did not converge");
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

    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux

#endif
