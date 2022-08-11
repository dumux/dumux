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
 * \ingroup PNMTwoPModel
 * \copydoc Dumux::PoreNetwork::TwoPNewtonSolver
 */
#ifndef DUMUX_PNM_NEWTON_SOLVER_HH
#define DUMUX_PNM_NEWTON_SOLVER_HH

#include <dumux/nonlinear/newtonsolver.hh>
#include "newtonconsistencychecks.hh"

namespace Dumux::PoreNetwork {
/*!
 * \ingroup PNMTwoPModel
 * \brief A two-phase PNM specific newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template<class Assembler, class LinearSolver, bool useRegularization = true>
class TwoPNewtonSolver;

template<class Assembler, class LinearSolver>
class TwoPNewtonSolver<Assembler, LinearSolver, true>: public Dumux::NewtonSolver<Assembler, LinearSolver>
{
    using ParentType =  Dumux::NewtonSolver<Assembler, LinearSolver>;
    using SolutionVector = typename Assembler::ResidualType;

public:
    using ParentType::ParentType;
    /*!
     * \brief Called after each Newton update
     *
     * \param uCurrentIter The current global solution vector
     * \param uLastIter The previous global solution vector
     */
     void newtonEnd(SolutionVector &uCurrentIter,
                    const SolutionVector &uLastIter) final
    {
        // call the method of the base class
        ParentType::newtonEnd(uCurrentIter, uLastIter);

        auto& gridVariables = this->assembler().gridVariables();
        auto& invasionState = gridVariables.gridFluxVarsCache().invasionState();

        // here the invasion state is updated
        invasionState.update(uCurrentIter, gridVariables.curGridVolVars(), gridVariables.gridFluxVarsCache());

        // If the solution is about to be accepted, check for accuracy and trigger a retry
        // with a decreased time step size if necessary.
        if (newtonConverged())
        {
            TwoPNewtonConsistencyChecks<typename Assembler::GridVariables, SolutionVector> checks;
            checks.performChecks(gridVariables, uCurrentIter, this->assembler().prevSol());
        }
    }

    /*!
     * \brief Returns true if the current solution can be considered to
     *        be accurate enough. We enforce an additional Newton iteration if the invasion
     *        switch was triggered in the last iteration.
     */
    bool newtonConverged() const final
    { return ParentType::newtonConverged(); }

    /*!
     * \brief Called if the Newton method broke down.
     * This method is called _after_ newtonEnd() and resets the invasion state.
     */
    void newtonFail(SolutionVector& u) final
    {
        ParentType::newtonFail(u);
        auto& gridVariables = this->assembler().gridVariables();
        gridVariables.gridFluxVarsCache().invasionState().reset();
    }

    /*!
     * \brief Called if the Newton method ended successfully
     * This method is called _after_ newtonEnd() and advances the invasion state.
     */
    void newtonSucceed() final
    {
        auto& gridVariables = this->assembler().gridVariables();
        gridVariables.gridFluxVarsCache().invasionState().advance();
    }
};
} // end namespace Dumux::PoreNetwork

#endif
