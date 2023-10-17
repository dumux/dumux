// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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

#include <dumux/linear/matrixconverter.hh>

namespace Dumux::Detail::LinearPDESolver {

template <class Solver, class Matrix>
using SetMatrixDetector = decltype(std::declval<Solver>().setMatrix(std::declval<std::shared_ptr<Matrix>>()));

template<class Solver, class Matrix>
static constexpr bool linearSolverHasSetMatrix()
{ return Dune::Std::is_detected<SetMatrixDetector, Solver, Matrix>::value; }

} // end namespace Dumux::Detail::LinearPDESolver

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
template <class Assembler, class LinearSolver,
          class Comm = Dune::Communication<Dune::MPIHelper::MPICommunicator>>
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
    using Communication = Comm;

    /*!
     * \brief The Constructor
     */
    LinearPDESolver(std::shared_ptr<Assembler> assembler,
                    std::shared_ptr<LinearSolver> linearSolver,
                    const Communication& comm = Dune::MPIHelper::getCommunication(),
                    const std::string& paramGroup = "")
    : ParentType(assembler, linearSolver)
    , paramGroup_(paramGroup)
    , reuseMatrix_(false)
    , comm_(comm)
    {
        initParams_(paramGroup);

        // set the linear system (matrix & residual) in the assembler
        this->assembler().setLinearSystem();
    }

    /*!
     * \brief The Constructor
     */
    LinearPDESolver(std::shared_ptr<Assembler> assembler,
                    std::shared_ptr<LinearSolver> linearSolver,
                    const std::string& paramGroup)
    : LinearPDESolver(assembler, linearSolver, Dune::MPIHelper::getCommunication(), paramGroup)
    {}

    /*!
     * \brief Solve a linear PDE system
     */
    bool apply(Variables& vars) override
    {
        Dune::Timer assembleTimer(false);
        Dune::Timer solveTimer(false);
        Dune::Timer updateTimer(false);

        if (verbose_ && enableDynamicOutput_)
            std::cout << "Assemble: r(x^k) = dS/dt + div F - q;   M = grad r"
                      << std::flush;

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

        if (verbose_ && enableDynamicOutput_)
            std::cout << "\rSolve: M deltax^k = r"
                      << clearRemainingLine << std::flush;

        // solve the resulting linear equation system
        solveTimer.start();

        // set the delta vector to zero
        ResidualVector deltaU = LinearAlgebraNativeBackend::zeros(Backend::size(Backend::dofs(vars)));

        // solve by calling the appropriate implementation depending on whether the linear solver
        // is capable of handling MultiType matrices or not
        bool converged = solveLinearSystem_(deltaU);
        solveTimer.stop();

        if (!converged)
            return false;

        ///////////////
        // update
        ///////////////
        if (verbose_ && enableDynamicOutput_)
            std::cout << "\rUpdate: x^(k+1) = x^k - deltax^k"
                      << clearRemainingLine << std::flush;

        // update the current solution and secondary variables
        updateTimer.start();
        auto uCurrent = Backend::dofs(vars);
        Backend::axpy(-1.0, deltaU, uCurrent);
        Backend::update(vars, uCurrent);
        if constexpr (!assemblerExportsVariables)
            this->assembler().updateGridVariables(Backend::dofs(vars));
        updateTimer.stop();

        if (verbose_)
        {
            const auto elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed();
            if (enableDynamicOutput_)
                std::cout << '\r';
            std::cout << "Assemble/solve/update time: "
                      <<  assembleTimer.elapsed() << "(" << 100*assembleTimer.elapsed()/elapsedTot << "%)/"
                      <<  solveTimer.elapsed() << "(" << 100*solveTimer.elapsed()/elapsedTot << "%)/"
                      <<  updateTimer.elapsed() << "(" << 100*updateTimer.elapsed()/elapsedTot << "%)"
                      << "\n";
        }

        return true;
    }

    /*!
     * \brief Solve a linear PDE system
     */
    void solve(Variables& vars) override
    {
        bool converged = apply(vars);
        if (!converged)
            DUNE_THROW(NumericalProblem, "Linear solver didn't converge.\n");
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
    {
        reuseMatrix_ = reuse;

        if constexpr (Detail::LinearPDESolver::linearSolverHasSetMatrix<LinearSolver, JacobianMatrix>())
            if (reuseMatrix_)
                this->linearSolver().setMatrix(this->assembler().jacobian());
    }

private:

    virtual bool solveLinearSystem_(ResidualVector& deltaU)
    {
        if constexpr (Detail::LinearPDESolver::linearSolverHasSetMatrix<LinearSolver, JacobianMatrix>())
        {
            if (reuseMatrix_)
                return this->linearSolver().solve(deltaU, this->assembler().residual());
        }
        else
        {
            if (reuseMatrix_ && comm_.size() > 1)
                DUNE_THROW(Dune::NotImplemented,
                    "Reuse matrix for parallel runs with a solver that doesn't support the setMatrix interface"
                );
        }

        assert(this->checkSizesOfSubMatrices(this->assembler().jacobian()) && "Matrix blocks have wrong sizes!");

        return this->linearSolver().solve(
            this->assembler().jacobian(),
            deltaU,
            this->assembler().residual()
        );
    }

    //! initialize the parameters by reading from the parameter tree
    void initParams_(const std::string& group = "")
    {
        verbose_ = (Dune::MPIHelper::getCommunication().rank() == 0);
        enableDynamicOutput_ = getParamFromGroup<bool>(group, "LinearPDESolver.EnableDynamicOutput", true);
    }

    //! switches on/off verbosity
    bool verbose_;

    //! further parameters
    bool enableDynamicOutput_;

    //! the parameter group for getting parameters from the parameter tree
    std::string paramGroup_;

    //! check if the matrix is supposed to be reused
    bool reuseMatrix_;

    //! The communication object (for distributed memory parallelism)
    Communication comm_;
};

} // end namespace Dumux

#endif
