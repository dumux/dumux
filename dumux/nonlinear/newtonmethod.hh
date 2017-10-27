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
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(Model);
NEW_PROP_TAG(NewtonController);
NEW_PROP_TAG(SolutionVector);
NEW_PROP_TAG(JacobianAssembler);
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
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) NewtonController;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;
public:
    NewtonMethod(Problem &problem)
        : problem_(problem)
    { }

    /*!
     * \brief Returns a reference to the current numeric problem.
     */
    Problem &problem()
    { return problem_; }

    /*!
     * \brief Returns a reference to the current numeric problem.
     */
    const Problem &problem() const
    { return problem_; }

    /*!
     * \brief Returns a reference to the numeric model.
     */
    Model &model()
    { return problem().model(); }

    /*!
     * \brief Returns a reference to the numeric model.
     */
    const Model &model() const
    { return problem().model(); }


    /*!
     * \brief Run the newton method. The controller is responsible
     *        for all the strategic decisions.
     */
    bool execute(NewtonController &ctl)
    {
        try {
            return execute_(ctl);
        }
        catch (const NumericalProblem &e) {
            if (ctl.verbose())
                std::cout << "Newton: Caught exception: \"" << e.what() << "\"\n";
            ctl.newtonFail();
            return false;
        }
    }

protected:
    bool execute_(NewtonController &ctl)
    {
        SolutionVector &uCurrentIter = model().curSol();
        SolutionVector &uLastIter = model().lastIter();
        SolutionVector deltaU(uCurrentIter);

        JacobianAssembler &jacobianAsm = model().jacobianAssembler();

        Dune::Timer assembleTimer(false);
        Dune::Timer solveTimer(false);
        Dune::Timer updateTimer(false);

        // tell the controller that we begin solving
        ctl.newtonBegin(*this, uCurrentIter);

        // execute the method as long as the controller thinks
        // that we should do another iteration
        while (ctl.newtonProceed(uCurrentIter))
        {
            // notify the controller that we're about to start
            // a new timestep
            ctl.newtonBeginStep();


            if (ctl.verbose()) {
                std::cout << "Assemble: r(x^k) = dS/dt + div F - q;   M = grad r";
                std::cout.flush();
            }

            ///////////////
            // assemble
            ///////////////

            // linearize the problem at the current solution
            assembleTimer.start();
            jacobianAsm.assemble();
            assembleTimer.stop();

            // make the current solution to the old one
            uLastIter = uCurrentIter;

            ///////////////
            // linear solve
            ///////////////

            // Clear the current line using an ansi escape
            // sequence.  for an explanation see
            // http://en.wikipedia.org/wiki/ANSI_escape_code
            const char clearRemainingLine[] = { 0x1b, '[', 'K', 0 };

            if (ctl.verbose()) {
                std::cout << "\rSolve: M deltax^k = r";
                std::cout << clearRemainingLine;
                std::cout.flush();
            }

            // solve the resulting linear equation system
            solveTimer.start();

            // set the delta vector to zero before solving the linear system!
            deltaU = 0;
            // ask the controller to solve the linearized system
            ctl.newtonSolveLinear(jacobianAsm.matrix(),
                                  deltaU,
                                  jacobianAsm.residual());
            solveTimer.stop();

            ///////////////
            // update
            ///////////////
            if (ctl.verbose()) {
                std::cout << "\rUpdate: x^(k+1) = x^k - deltax^k";
                std::cout << clearRemainingLine;
                std::cout.flush();
            }

            updateTimer.start();
            // update the current solution (i.e. uOld) with the delta
            // (i.e. u). The result is stored in u
            ctl.newtonUpdate(uCurrentIter, uLastIter, deltaU);
            updateTimer.stop();

            // tell the controller that we're done with this iteration
            ctl.newtonEndStep(uCurrentIter, uLastIter);
        }

        // tell the controller that we're done
        ctl.newtonEnd();

        if (ctl.verbose()) {
            Scalar elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed();
            std::cout << "Assemble/solve/update time: "
                      <<  assembleTimer.elapsed() << "(" << 100*assembleTimer.elapsed()/elapsedTot << "%)/"
                      <<  solveTimer.elapsed() << "(" << 100*solveTimer.elapsed()/elapsedTot << "%)/"
                      <<  updateTimer.elapsed() << "(" << 100*updateTimer.elapsed()/elapsedTot << "%)"
                      << "\n";
        }

        if (!ctl.newtonConverged()) {
            ctl.newtonFail();
            return false;
        }

        ctl.newtonSucceed();
        return true;
    }

private:
    Problem &problem_;
};

}

#endif
