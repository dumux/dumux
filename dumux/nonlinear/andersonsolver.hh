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
 * \ingroup Common
 * \brief Defines a high-level interface for a PDESolver
 */
#ifndef DUMUX_NONLINEAR_ANDERSONSOLVER_HH
#define DUMUX_NONLINEAR_ANDERSONSOLVER_HH

#include <memory>

#include <dumux/common/timeloop.hh>
#include <dumux/common/pdesolver.hh>

namespace Dumux {

template<class Assembler, class LinearSolver>
class AndersonSolver: public PDESolver<Assembler, LinearSolver>
{
    using ParentType = PDESolver<Assembler, LinearSolver>;
    using SolutionVector = typename Assembler::ResidualType;
    using Scalar = typename Assembler::Scalar;
    using TimeLoop = TimeLoopBase<Scalar>;

public:
    AndersonSolver(std::shared_ptr<Assembler> assembler,
                   std::shared_ptr<LinearSolver> linearSolver,
                   const std::string& paramGroup = "")
    : ParentType(assembler, linearSolver)
    {
        initParams_(paramGroup);

        this->assembler().setLinearSystem();
    }


    /*!
     * \brief Solve the given PDE system (usually assemble + solve linear system + update)
     * \param sol a solution vector possbilty containing an initial solution
     */
    void solve(SolutionVector& sol)
    {}

    void solve(SolutionVector& uCurrentIter, TimeLoop& timeLoop) override
    {
        if (this->assembler().isStationaryProblem())
            DUNE_THROW(Dune::InvalidStateException, "Using time step control with stationary problem makes no sense!");

        // try solving the non-linear system
        for (std::size_t i = 0; i <= maxTimeStepDivisions_; ++i)
        {
            // linearize & solve
            const bool converged = solve_(uCurrentIter);

            if (converged)
                return;

            else if (!converged && i < maxTimeStepDivisions_)
            {
                // set solution to previous solution
                uCurrentIter = this->assembler().prevSol();

                // reset the grid variables to the previous solution
                this->assembler().resetTimeStep(uCurrentIter);

                if (verbosity_ >= 1)
                    std::cout << "Anderson solver did not converge with dt = "
                              << timeLoop.timeStepSize() << " seconds. Retrying with time step of "
                              << timeLoop.timeStepSize() * retryTimeStepReductionFactor_ << " seconds\n";

                // try again with dt = dt * retryTimeStepReductionFactor_
                timeLoop.setTimeStepSize(timeLoop.timeStepSize() * retryTimeStepReductionFactor_);
            }

            else
            {
                DUNE_THROW(NumericalProblem, "Anderson solver didn't converge after "
                                             << maxTimeStepDivisions_ << " time-step divisions. dt="
                                             << timeLoop.timeStepSize() << '\n');
            }
        }
    }

private:
    bool solve_(SolutionVector& uCurrentIter)
    {
        // the given solution is the initial guess
        SolutionVector uLastIter(uCurrentIter);
        SolutionVector deltaU(uCurrentIter);

        Dune::Timer assembleTimer(false);
        Dune::Timer solveTimer(false);
        Dune::Timer updateTimer(false);

        try
        {
            newtonBegin(uCurrentIter);

            // execute the method as long as the solver thinks
            // that we should do another iteration
            bool converged = false;
            while (newtonProceed(uCurrentIter, converged))
            {
                // notify the solver that we're about to start
                // a new timestep
                newtonBeginStep(uCurrentIter);

                // make the current solution to the old one
                if (numSteps_ > 0)
                    uLastIter = uCurrentIter;

                if (verbosity_ >= 1 && enableDynamicOutput_)
                    std::cout << "Assemble: r(x^k) = dS/dt + div F - q;   M = grad r"
                              << std::flush;

                ///////////////
                // assemble
                ///////////////

                // linearize the problem at the current solution
                assembleTimer.start();
                assembleLinearSystem(uCurrentIter);
                assembleTimer.stop();

                ///////////////
                // linear solve
                ///////////////

                // Clear the current line using an ansi escape
                // sequence.  for an explanation see
                // http://en.wikipedia.org/wiki/ANSI_escape_code
                const char clearRemainingLine[] = { 0x1b, '[', 'K', 0 };

                if (verbosity_ >= 1 && enableDynamicOutput_)
                    std::cout << "\rSolve: M deltax^k = r"
                              << clearRemainingLine << std::flush;

                // solve the resulting linear equation system
                solveTimer.start();

                // set the delta vector to zero before solving the linear system!
                deltaU = 0;

                solveLinearSystem(deltaU);
                solveTimer.stop();

                ///////////////
                // update
                ///////////////
                if (verbosity_ >= 1 && enableDynamicOutput_)
                    std::cout << "\rUpdate: x^(k+1) = x^k - deltax^k"
                              << clearRemainingLine << std::flush;

                updateTimer.start();
                // update the current solution (i.e. uOld) with the delta
                // (i.e. u). The result is stored in u
                newtonUpdate(uCurrentIter, uLastIter, deltaU);
                updateTimer.stop();

                // tell the solver that we're done with this iteration
                newtonEndStep(uCurrentIter, uLastIter);

                // if a convergence writer was specified compute residual and write output
                if (convergenceWriter_)
                {
                    this->assembler().assembleResidual(uCurrentIter);
                    convergenceWriter_->write(uCurrentIter, deltaU, this->assembler().residual());
                }

                // detect if the method has converged
                converged = newtonConverged();
            }

            // tell solver we are done
            newtonEnd();

            // reset state if Newton failed
            if (!newtonConverged())
            {
                totalWastedIter_ += numSteps_;
                newtonFail(uCurrentIter);
                return false;
            }

            totalSucceededIter_ += numSteps_;
            numConverged_++;

            // tell solver we converged successfully
            newtonSucceed();

            if (verbosity_ >= 1) {
                const auto elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed();
                std::cout << "Assemble/solve/update time: "
                          <<  assembleTimer.elapsed() << "(" << 100*assembleTimer.elapsed()/elapsedTot << "%)/"
                          <<  solveTimer.elapsed() << "(" << 100*solveTimer.elapsed()/elapsedTot << "%)/"
                          <<  updateTimer.elapsed() << "(" << 100*updateTimer.elapsed()/elapsedTot << "%)"
                          << "\n";
            }
            return true;

        }
        catch (const NumericalProblem &e)
        {
            if (verbosity_ >= 1)
                std::cout << "Newton: Caught exception: \"" << e.what() << "\"\n";

            totalWastedIter_ += numSteps_;
            newtonFail(uCurrentIter);
            return false;
        }
    }

    //! initialize the parameters by reading from the parameter tree
    void initParams_(const std::string& group = "")
    {
        enableAbsoluteResidualCriterion_ = getParamFromGroup<bool>(group, "Anderson.EnableAbsoluteResidualCriterion");
        enableShiftCriterion_ = getParamFromGroup<bool>(group, "Anderson.EnableShiftCriterion");
        enableResidualCriterion_ = getParamFromGroup<bool>(group, "Anderson.EnableResidualCriterion") || enableAbsoluteResidualCriterion_;
        satisfyResidualAndShiftCriterion_ = getParamFromGroup<bool>(group, "Anderson.SatisfyResidualAndShiftCriterion");

        setMaxRelativeShift(getParamFromGroup<Scalar>(group, "Anderson.MaxRelativeShift"));
        setMaxAbsoluteResidual(getParamFromGroup<Scalar>(group, "Anderson.MaxAbsoluteResidual"));
        setResidualReduction(getParamFromGroup<Scalar>(group, "Anderson.ResidualReduction"));
        setTargetSteps(getParamFromGroup<int>(group, "Anderson.TargetSteps"));
        setMinSteps(getParamFromGroup<int>(group, "Anderson.MinSteps"));
        setMaxSteps(getParamFromGroup<int>(group, "Anderson.MaxSteps"));

        maxTimeStepDivisions_ = getParamFromGroup<std::size_t>(group, "Anderson.MaxTimeStepDivisions", 10);
        retryTimeStepReductionFactor_ = getParamFromGroup<Scalar>(group, "Anderson.RetryTimeStepReductionFactor", 0.5);

        verbosity_ = getParamFromGroup<int>(group, "Anderson.Verbosity", 2);
        numSteps_ = 0;
    }

    //! optimal number of iterations we want to achieve
    int targetSteps_;
    //! minimum number of iterations we do
    int minSteps_;
    //! maximum number of iterations we do before giving up
    int maxSteps_;
    //! actual number of steps done so far
    int numSteps_;

    // residual criterion variables
    Scalar reduction_;
    Scalar residualNorm_;
    Scalar lastReduction_;
    Scalar initialResidual_;

    // shift criterion variables
    Scalar shift_;
    Scalar lastShift_;

    //! message stream to be displayed at the end of iterations
    std::ostringstream endIterMsgStream_;
    //! the verbosity level
    int verbosity_;

    Scalar shiftTolerance_;
    Scalar reductionTolerance_;
    Scalar residualTolerance_;

    // time step control
    std::size_t maxTimeStepDivisions_;
    Scalar retryTimeStepReductionFactor_;

    // further parameters
    bool useLineSearch_;
    bool useChop_;
    bool enableAbsoluteResidualCriterion_;
    bool enableShiftCriterion_;
    bool enableResidualCriterion_;
    bool satisfyResidualAndShiftCriterion_;
    bool enableDynamicOutput_;

    //! the parameter group for getting parameters from the parameter tree
    std::string paramGroup_;
};

} // namespace Dumux

#endif
