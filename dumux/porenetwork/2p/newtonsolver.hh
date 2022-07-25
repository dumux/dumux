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
template<class Assembler, class LinearSolver,
         template<class, class> class NewtonConsistencyChecks = TwoPNewtonConsistencyChecks>
class TwoPNewtonSolver : public Dumux::NewtonSolver<Assembler, LinearSolver>
{
    using ParentType =  Dumux::NewtonSolver<Assembler, LinearSolver>;
    using ConstSolutionVector = typename Assembler::ResidualType;
    using SolutionVector = std::remove_const_t<ConstSolutionVector>;
    using Variables = typename PDESolver<Assembler, LinearSolver>::Variables;
    using Backend = VariablesBackend<Variables>;

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

        invasionState.update(uCurrentIter, gridVariables.curGridVolVars(), gridVariables.gridFluxVarsCache());
    }

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

    // /*!
    //  * \brief Called if the Newton method ended successfully
    //  * This method is called _after_ newtonEnd() and advances the invasion state.
    //  */
    void newtonSucceed() final
    {
        auto& gridVariables = this->assembler().gridVariables();
        gridVariables.gridFluxVarsCache().invasionState().advance();
    }

    void setIndex(std::vector<bool> poreindex)
    { preIndex_ = poreindex; }

private:
    std::vector<bool> preIndex_; // return the pore indicies which are predicted to be invaded
};

} // end namespace Dumux::PoreNetwork

#endif
