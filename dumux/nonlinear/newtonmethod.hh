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
 * \brief The algorithmic part of the multi dimensional newton method.
 *
 * In order to use the method you need a Newtoncontroller
 */
#ifndef DUMUX_NEWTONMETHOD_HH
#define DUMUX_NEWTONMETHOD_HH

#include <memory>
#include <iostream>

#include <dune/common/timer.hh>
#include <dune/istl/istlexception.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dune/common/deprecated.hh>

#warning "This file is deprecated. Use NewtonSolver instead."

namespace Dumux {

/*!
 * \ingroup Newton
 * \ingroup Nonlinear
 * \brief The algorithmic part of the multi dimensional newton method.
 *
 * \tparam NewtonController the controller class is the driver implementation
 * \tparam JacobianAssembler an assembler for the Jacobian and the residual
 * \tparam LinearSolver the linear solver used to solve one iteration
 */
template <class NewtonController, class JacobianAssembler, class LinearSolver>
class DUNE_DEPRECATED_MSG("Use NewtonSolver instead.") NewtonMethod
{
    //! provide an interface as a form of type erasure
    //! this is the minimal requirements a convergence write passed to a newton method has to fulfill
    struct ConvergenceWriterInferface
    {
        template<class SolutionVector>
        void write(const SolutionVector &uLastIter, const SolutionVector &deltaU, const SolutionVector &residual) {}
    };

public:
    NewtonMethod(std::shared_ptr<NewtonController> controller,
                 std::shared_ptr<JacobianAssembler> assembler,
                 std::shared_ptr<LinearSolver> linearSolver,
                 const std::string& modelParamGroup = "")
    : controller_(controller)
    , assembler_(assembler)
    , linearSolver_(linearSolver)
    {
        // set the linear system (matrix & residual) in the assembler
        assembler_->setLinearSystem();

        // set a different default for the linear solver residual reduction
        // within the Newton the linear solver doesn't need to solve too exact
        using Scalar = typename LinearSolver::Scalar;
        linearSolver_->setResidualReduction(getParamFromGroup<Scalar>(modelParamGroup, "LinearSolver.ResidualReduction", 1e-6));
    }

    /*!
     * \brief Run the newton method to solve a non-linear system.
     *        The controller is responsible for all the strategic decisions.
     */
    template<class SolutionVector, class ConvergenceWriter = ConvergenceWriterInferface>
    bool solve(SolutionVector& uCurrentIter, const std::unique_ptr<ConvergenceWriter>& convWriter = nullptr)
    {
        try
        {
            // the given solution is the initial guess
            SolutionVector uLastIter(uCurrentIter);
            SolutionVector deltaU(uCurrentIter);

            Dune::Timer assembleTimer(false);
            Dune::Timer solveTimer(false);
            Dune::Timer updateTimer(false);

            // tell the controller that we begin solving
            controller_->newtonBegin(uCurrentIter);

            // execute the method as long as the controller thinks
            // that we should do another iteration
            while (controller_->newtonProceed(uCurrentIter, controller_->newtonConverged()))
            {
                // notify the controller that we're about to start
                // a new timestep
                controller_->newtonBeginStep(uCurrentIter);

                // make the current solution to the old one
                if (controller_->newtonNumSteps() > 0)
                    uLastIter = uCurrentIter;

                if (controller_->verbose()) {
                    std::cout << "Assemble: r(x^k) = dS/dt + div F - q;   M = grad r";
                    std::cout.flush();
                }

                ///////////////
                // assemble
                ///////////////

                // linearize the problem at the current solution
                assembleTimer.start();
                controller_->assembleLinearSystem(*assembler_, uCurrentIter);
                assembleTimer.stop();

                ///////////////
                // linear solve
                ///////////////

                // Clear the current line using an ansi escape
                // sequence.  for an explanation see
                // http://en.wikipedia.org/wiki/ANSI_escape_code
                const char clearRemainingLine[] = { 0x1b, '[', 'K', 0 };

                if (controller_->verbose()) {
                    std::cout << "\rSolve: M deltax^k = r";
                    std::cout << clearRemainingLine;
                    std::cout.flush();
                }

                // solve the resulting linear equation system
                solveTimer.start();

                // set the delta vector to zero before solving the linear system!
                deltaU = 0;
                // ask the controller to solve the linearized system
                controller_->solveLinearSystem(*linearSolver_,
                                               assembler_->jacobian(),
                                               deltaU,
                                               assembler_->residual());
                solveTimer.stop();

                ///////////////
                // update
                ///////////////
                if (controller_->verbose()) {
                    std::cout << "\rUpdate: x^(k+1) = x^k - deltax^k";
                    std::cout << clearRemainingLine;
                    std::cout.flush();
                }

                updateTimer.start();
                // update the current solution (i.e. uOld) with the delta
                // (i.e. u). The result is stored in u
                controller_->newtonUpdate(*assembler_, uCurrentIter, uLastIter, deltaU);
                updateTimer.stop();

                // tell the controller that we're done with this iteration
                controller_->newtonEndStep(*assembler_, uCurrentIter, uLastIter);

                // if a convergence writer was specified compute residual and write output
                if (convWriter)
                {
                    assembler_->assembleResidual(uCurrentIter);
                    convWriter->write(uLastIter, deltaU, assembler_->residual());
                }
            }

            // tell controller we are done
            controller_->newtonEnd();

            // reset state if newton failed
            if (!controller_->newtonConverged())
            {
                controller_->newtonFail(*assembler_, uCurrentIter);
                return false;
            }

            // tell controller we converged successfully
            controller_->newtonSucceed();

            if (controller_->verbose()) {
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
            if (controller_->verbose())
                std::cout << "Newton: Caught exception: \"" << e.what() << "\"\n";
            controller_->newtonFail(*assembler_, uCurrentIter);
            return false;
        }
    }

private:
    std::shared_ptr<NewtonController> controller_;
    std::shared_ptr<JacobianAssembler> assembler_;
    std::shared_ptr<LinearSolver> linearSolver_;
};

} // end namespace Dumux

#endif
