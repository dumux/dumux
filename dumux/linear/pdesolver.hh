// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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
 * \ingroup Linear
 * \brief A high-level solver interface for a linear PDE solver
 */
#ifndef DUMUX_LINEAR_PDE_SOLVER_HH
#define DUMUX_LINEAR_PDE_SOLVER_HH

#include <cmath>
#include <memory>
#include <iostream>
#include <type_traits>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpicommunication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/typetraits/vector.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/pdesolver.hh>
#include <dumux/common/variablesbackend.hh>

// remove after deprecated code is removed (after 3.7)
#define DUMUX_SUPPRESS_LINEAR_SOLVER_ACCEPTS_MULTITYPEMATRIX_WARNING
#include <dumux/linear/linearsolveracceptsmultitypematrix.hh>
#undef DUMUX_SUPPRESS_LINEAR_SOLVER_ACCEPTS_MULTITYPEMATRIX_WARNING

#include <dumux/linear/matrixconverter.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief An implementation of a linear PDE solver
 * \tparam Assembler the assembler
 * \tparam LinearSolver the linear solver
 * \tparam Comm the communication object used to communicate with all processes
 * \note If you want to specialize only some methods but are happy with the
 *       defaults of the reference solver, derive your solver from
 *       this class and simply overload the required methods.
 */
template <class Assembler, class LinearSolver>
class LinearPDESolver : public PDESolver<Assembler, LinearSolver>
{
    using ParentType = PDESolver<Assembler, LinearSolver>;
    using Scalar = typename Assembler::Scalar;
    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using SolutionVector = typename Assembler::SolutionVector;
    using ResidualVector = typename Assembler::ResidualType;
    using TimeLoop = TimeLoopBase<Scalar>;
    using Backend = VariablesBackend<typename ParentType::Variables>;
    using LinearAlgebraNativeBackend = VariablesBackend<ResidualVector>;
    static constexpr bool assemblerExportsVariables = Detail::PDESolver::assemblerExportsVariables<Assembler>;

public:
    using typename ParentType::Variables;

    /*!
     * \brief The Constructor
     */
    LinearPDESolver(std::shared_ptr<Assembler> assembler,
                    std::shared_ptr<LinearSolver> linearSolver,
                    const std::string& paramGroup = "")
    : ParentType(assembler, linearSolver)
    , paramGroup_(paramGroup)
    , reuseMatrix_(false)
    {
        initParams_(paramGroup);

        // set the linear system (matrix & residual) in the assembler
        this->assembler().setLinearSystem();
    }

    /*!
     * \brief Solve a linear PDE system
     */
    void solve(Variables& vars) override
    {
        Dune::Timer assembleTimer(false);
        Dune::Timer solveTimer(false);
        Dune::Timer updateTimer(false);

        if (verbose_) {
            std::cout << "Assemble: r(x^k) = dS/dt + div F - q;   M = grad r"
                      << std::flush;
        }

        ///////////////
        // assemble
        ///////////////

        // linearize the problem at the current solution
        assembleTimer.start();
        if (reuseMatrix_)
            this->assembler().assembleResidual(vars);
        else
            this->assembler().assembleJacobianAndResidual(vars);
        assembleTimer.stop();

        ///////////////
        // linear solve
        ///////////////

        // Clear the current line using an ansi escape
        // sequence.  for an explanation see
        // http://en.wikipedia.org/wiki/ANSI_escape_code
        const char clearRemainingLine[] = { 0x1b, '[', 'K', 0 };

        if (verbose_) {
            std::cout << "\rSolve: M deltax^k = r";
            std::cout << clearRemainingLine
                      << std::flush;
        }

        // solve the resulting linear equation system
        solveTimer.start();

        // set the delta vector to zero
        ResidualVector deltaU = LinearAlgebraNativeBackend::zeros(Backend::size(Backend::dofs(vars)));

        // solve by calling the appropriate implementation depending on whether the linear solver
        // is capable of handling MultiType matrices or not
        bool converged = solveLinearSystem_(deltaU);
        solveTimer.stop();

        if (!converged)
            DUNE_THROW(NumericalProblem, "Linear solver didn't converge.\n");

        ///////////////
        // update
        ///////////////
        if (verbose_) {
            std::cout << "\rUpdate: x^(k+1) = x^k - deltax^k";
            std::cout << clearRemainingLine;
            std::cout.flush();
        }

        // update the current solution and secondary variables
        updateTimer.start();
        auto uCurrent = Backend::dofs(vars);
        Backend::axpy(-1.0, deltaU, uCurrent);
        Backend::update(vars, uCurrent);
        if constexpr (!assemblerExportsVariables)
            this->assembler().updateGridVariables(Backend::dofs(vars));
        updateTimer.stop();

        if (verbose_) {
            const auto elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed();
            std::cout << "Assemble/solve/update time: "
                      <<  assembleTimer.elapsed() << "(" << 100*assembleTimer.elapsed()/elapsedTot << "%)/"
                      <<  solveTimer.elapsed() << "(" << 100*solveTimer.elapsed()/elapsedTot << "%)/"
                      <<  updateTimer.elapsed() << "(" << 100*updateTimer.elapsed()/elapsedTot << "%)"
                      << "\n";
        }
    }

    /*!
     * \brief output statistics / report
     * \todo Implement some solver statistics output
     */
    void report(std::ostream& sout = std::cout) const
    {}

    /*!
     * \brief Suggest a new time-step size based on the old time-step size.
     * \note For compatibility with other PDE solvers (e.g. Newton)
     */
    Scalar suggestTimeStepSize(Scalar oldTimeStep) const
    {
        return oldTimeStep;
    }

    /*!
     * \brief Specifies if the solver ought to be chatty.
     */
    void setVerbose(bool val)
    { verbose_ = val; }

    /*!
     * \brief Returns true if the solver ought to be chatty.
     */
    bool verbose() const
    { return verbose_ ; }

    /*!
     * \brief Returns the parameter group
     */
    const std::string& paramGroup() const
    { return paramGroup_; }

    /*!
     * \brief Set whether the matrix should be reused
     * \note If this is set to true, the matrix will not be assembled. Make
     *       sure there is an assembled matrix that can be reused before
     *       setting this flag to true.
     */
    void reuseMatrix(bool reuse = true)
    { reuseMatrix_ = reuse; }

private:

    virtual bool solveLinearSystem_(ResidualVector& deltaU)
    {
        return solveLinearSystemImpl_(this->linearSolver(),
                                      this->assembler().jacobian(),
                                      deltaU,
                                      this->assembler().residual());
    }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * Specialization for linear solvers that can handle MultiType matrices.
     *
     */
    template<class V = ResidualVector>
    typename std::enable_if_t<!isMultiTypeBlockVector<V>(), bool>
    solveLinearSystemImpl_(LinearSolver& ls,
                           JacobianMatrix& A,
                           ResidualVector& x,
                           ResidualVector& b)
    {
        return ls.solve(A, x, b);
    }


    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * Specialization for linear solvers that can handle MultiType matrices.
     *
     */
    template<class LS = LinearSolver, class V = ResidualVector>
    typename std::enable_if_t<linearSolverAcceptsMultiTypeMatrix<LS>() &&
                              isMultiTypeBlockVector<V>(), bool>
    solveLinearSystemImpl_(LinearSolver& ls,
                           JacobianMatrix& A,
                           ResidualVector& x,
                           ResidualVector& b)
    {
        assert(this->checkSizesOfSubMatrices(A) && "Sub-blocks of MultiTypeBlockMatrix have wrong sizes!");
        return ls.solve(A, x, b);
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
    template<class LS = LinearSolver, class V = ResidualVector>
    [[deprecated("After 3.7 LinearPDESolver will no longer support conversion of multitype matrices for solvers that don't support this feature!")]]
    typename std::enable_if_t<!linearSolverAcceptsMultiTypeMatrix<LS>() &&
                              isMultiTypeBlockVector<V>(), bool>
    solveLinearSystemImpl_(LinearSolver& ls,
                           JacobianMatrix& A,
                           ResidualVector& x,
                           ResidualVector& b)
    {
        assert(this->checkSizesOfSubMatrices(A) && "Sub-blocks of MultiTypeBlockMatrix have wrong sizes!");

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

    //! initialize the parameters by reading from the parameter tree
    void initParams_(const std::string& group = "")
    {
        verbose_ = (Dune::MPIHelper::getCommunication().rank() == 0);
    }

    //! switches on/off verbosity
    bool verbose_;

    //! the parameter group for getting parameters from the parameter tree
    std::string paramGroup_;

    //! check if the matrix is supposed to be reused
    bool reuseMatrix_;
};

} // end namespace Dumux

#endif
