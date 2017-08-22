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
 *
 * \brief The algorithmic part of the multi dimensional newton method.
 *
 * In order to use the method you need a NewtonController.
 */
#ifndef DUMUX_NEWTONMETHOD_HH
#define DUMUX_NEWTONMETHOD_HH

#include <dumux/common/exceptions.hh>
#include <dumux/common/propertysystem.hh>

#include <dune/common/timer.hh>
#include <dune/istl/istlexception.hh>

#include <iostream>

namespace Dumux
{
// forward declaration of property tags
namespace Properties
{
// create a new type tag for models which apply the newton method
NEW_TYPE_TAG(NewtonMethod);

NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(NewtonController);
NEW_PROP_TAG(SolutionVector);
NEW_PROP_TAG(JacobianAssembler);
NEW_PROP_TAG(JacobianMatrix);
NEW_PROP_TAG(WriteConvergence);
}

/*!
 * \ingroup Newton
 * \brief The algorithmic part of the multi dimensional newton method.
 *
 * In order to use the method you need a NewtonController.
 */
template <class TypeTag>
class NewtonMethod
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

public:
    NewtonMethod()
    : matrix_(std::make_shared<JacobianMatrix>())
    , residual_(std::make_shared<SolutionVector>())
    {}

    /*!
     * \brief Run the newton method to solve a non-linear system.
     *        The controller is responsible for all the strategic decisions.
     */
    template<class NewtonController, class JacobianAssembler, class LinearSolver, class SolutionVector>
    bool solve(NewtonController& controller, JacobianAssembler& assembler, LinearSolver& linearSolver,
               SolutionVector& u, const SolutionVector& uPrev)
    {
        try
        {
            // set the linear system (matrix & residual) in the assembler
            assembler.setLinearSystem(matrix_, residual_);

            // the current solution is the initial guess
            SolutionVector& uCurrentIter = u;
            SolutionVector uLastIter(uCurrentIter);
            SolutionVector deltaU(uCurrentIter);

            Dune::Timer assembleTimer(false);
            Dune::Timer solveTimer(false);
            Dune::Timer updateTimer(false);

            // tell the controller that we begin solving
            controller.newtonBegin(uCurrentIter);

            // execute the method as long as the controller thinks
            // that we should do another iteration
            while (controller.newtonProceed(uCurrentIter))
            {
                // notify the controller that we're about to start
                // a new timestep
                controller.newtonBeginStep();

                // make the current solution to the old one
                if (controller.newtonNumSteps() > 0)
                    uLastIter = uCurrentIter;

                if (controller.verbose()) {
                    std::cout << "Assemble: r(x^k) = dS/dt + div F - q;   M = grad r";
                    std::cout.flush();
                }

                ///////////////
                // assemble
                ///////////////

                // linearize the problem at the current solution
                assembleTimer.start();
                controller.assembleLinearSystem(assembler, u, uPrev);
                assembleTimer.stop();

                ///////////////
                // linear solve
                ///////////////

                // Clear the current line using an ansi escape
                // sequence.  for an explanation see
                // http://en.wikipedia.org/wiki/ANSI_escape_code
                const char clearRemainingLine[] = { 0x1b, '[', 'K', 0 };

                if (controller.verbose()) {
                    std::cout << "\rSolve: M deltax^k = r";
                    std::cout << clearRemainingLine;
                    std::cout.flush();
                }

                // solve the resulting linear equation system
                solveTimer.start();

                // set the delta vector to zero before solving the linear system!
                deltaU = 0;
                // ask the controller to solve the linearized system
                controller.solveLinearSystem(linearSolver,
                                             matrix(),
                                             deltaU,
                                             residual());
                solveTimer.stop();

                ///////////////
                // update
                ///////////////
                if (controller.verbose()) {
                    std::cout << "\rUpdate: x^(k+1) = x^k - deltax^k";
                    std::cout << clearRemainingLine;
                    std::cout.flush();
                }

                updateTimer.start();
                // update the current solution (i.e. uOld) with the delta
                // (i.e. u). The result is stored in u
                controller.newtonUpdate(assembler, uCurrentIter, uLastIter, deltaU, uPrev);
                updateTimer.stop();

                // tell the controller that we're done with this iteration
                controller.newtonEndStep(assembler, uCurrentIter, uLastIter);
            }

            // tell the controller that we're done
            controller.newtonEnd();

            if (controller.verbose()) {
                Scalar elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed();
                std::cout << "Assemble/solve/update time: "
                          <<  assembleTimer.elapsed() << "(" << 100*assembleTimer.elapsed()/elapsedTot << "%)/"
                          <<  solveTimer.elapsed() << "(" << 100*solveTimer.elapsed()/elapsedTot << "%)/"
                          <<  updateTimer.elapsed() << "(" << 100*updateTimer.elapsed()/elapsedTot << "%)"
                          << "\n";
            }

            if (!controller.newtonConverged())
            {
                controller.newtonFail();
                return false;
            }

            controller.newtonSucceed();
            return true;

        }
        catch (const NumericalProblem &e)
        {
            if (controller.verbose())
                std::cout << "Newton: Caught exception: \"" << e.what() << "\"\n";
            controller.newtonFail();
            return false;
        }
    }

    //! The jacobian for the Newton method
    JacobianMatrix& matrix()
    { return *matrix_; }

    //! The residual for the Newton method
    SolutionVector& residual()
    { return *residual_; }

private:

    std::shared_ptr<JacobianMatrix> matrix_; //! The jacobian for the Newton method
    std::shared_ptr<SolutionVector> residual_; //! The residual for the Newton method
};

} // end namespace Dumux

#endif
