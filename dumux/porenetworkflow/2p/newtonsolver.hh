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
 * \copydoc Dumux::PNMNewtonSolver
 */
#ifndef DUMUX_PNM_NEWTON_SOLVER_HH
#define DUMUX_PNM_NEWTON_SOLVER_HH

#include <dumux/nonlinear/newtonsolver.hh>
#include "newtonconsistencychecks.hh"

namespace Dumux {
/*!
 * \ingroup PNMTwoP
 * \brief A two-phase PNM specific  newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template<class Assembler, class LinearSolver,
         template<class, class> class NewtonConsistencyChecks = PNMNewtonConsistencyChecks>
class PNMNewtonSolver : public NewtonSolver<Assembler, LinearSolver>
{
    using ParentType =  NewtonSolver<Assembler, LinearSolver>;
    using SolutionVector = typename Assembler::ResidualType;

public:
    using ParentType::ParentType;

    /*!
     * \brief Called after each Newton update
     *
     * \param uCurrentIter The current global solution vector
     * \param uLastIter The previous global solution vector
     */
     void newtonEndStep(SolutionVector &uCurrentIter,
                        const SolutionVector &uLastIter) final
    {
        // call the method of the base class
        ParentType::newtonEndStep(uCurrentIter, uLastIter);

        auto& gridVariables = this->assembler().gridVariables();
        auto& invasionState = gridVariables.gridFluxVarsCache().invasionState();
        switchedInLastIteration_ = invasionState.update(uCurrentIter, gridVariables.curGridVolVars(), gridVariables.gridFluxVarsCache());

        // If the solution is about to be accepted, check for accuracy and trigger a retry
        // with a decreased time step size if necessary.
        if (newtonConverged())
        {
            NewtonConsistencyChecks<typename Assembler::GridVariables, SolutionVector> checks;
            checks.performChecks(gridVariables, uCurrentIter, this->assembler().prevSol());
        }
    }

    /*!
     * \brief Returns true if the current solution can be considered to
     *        be accurate enough
     */
    bool newtonConverged() const final
    {
        if (switchedInLastIteration_)
            return false;

        return ParentType::newtonConverged();
    }

    void newtonFail(SolutionVector& u) final
    {
        ParentType::newtonFail(u);
        auto& gridVariables = this->assembler().gridVariables();
        gridVariables.gridFluxVarsCache().invasionState().reset();
    }

    /*!
     * \brief Called if the Newton method ended successfully
     * This method is called _after_ newtonEnd()
     */
    void newtonSucceed() final
    {
        auto& gridVariables = this->assembler().gridVariables();
        gridVariables.gridFluxVarsCache().invasionState().advance();
    }

private:
    bool switchedInLastIteration_{false};
};
}

#endif
