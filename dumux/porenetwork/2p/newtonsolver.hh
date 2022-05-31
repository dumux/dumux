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

        auto& gridVariables = this->assemblerReg().gridVariables();
        auto& invasionState = gridVariables.gridFluxVarsCache().invasionState();

        invasionState.update(uCurrentIter, gridVariables.curGridVolVars(), gridVariables.gridFluxVarsCache());

        // If the solution is about to be accepted, check for accuracy and trigger a retry
        // with a decreased time step size if necessary.
        if (newtonConverged())
        {
            NewtonConsistencyChecks<typename Assembler::GridVariables, SolutionVector> checks;
            checks.performChecks(gridVariables, uCurrentIter, this->assemblerReg().prevSol());
        }
    }

    /*!
     * \brief Returns true if the current solution can be considered to
     *        be accurate enough. We enforce an additional Newton iteration if the invasion
     *        switch was triggered in the last iteration.
     */
    bool newtonConverged() const final
    {
        //!TODO: last iteration not chopped
        if (lastIterationWasChopped_)
            return false;

        return ParentType::newtonConverged();
    }

    /*!
     * \brief Called if the Newton method broke down.
     * This method is called _after_ newtonEnd() and resets the invasion state.
     */
    void newtonFail(SolutionVector& u) final
    {
        ParentType::newtonFail(u);
        auto& gridVariables = this->assemblerReg().gridVariables();
        gridVariables.gridFluxVarsCache().invasionState().reset();
        lastIterationWasChopped_ = false;
    }

    // /*!
    //  * \brief Called if the Newton method ended successfully
    //  * This method is called _after_ newtonEnd() and advances the invasion state.
    //  */
    void newtonSucceed() final
    {
        auto& gridVariables = this->assemblerReg().gridVariables();
        gridVariables.gridFluxVarsCache().invasionState().advance();
    }

    void setIndex(std::vector<bool> poreindex)
    { preIndex_ = poreindex; }

private:

    /*!
     * \brief Update the current solution of the Newton method
     *
     * \param varsCurrentIter The variables after the current Newton iteration \f$ u^{k+1} \f$
     * \param uLastIter The solution after the last Newton iteration \f$ u^k \f$
     * \param deltaU The vector of differences between the last
     *               iterative solution and the next one \f$ \Delta u^k \f$
     */
    void choppedUpdate_(Variables &varsCurrentIter,
                        const SolutionVector &uLastIter,
                        const SolutionVector &deltaU) final
    {
        auto uCurrentIter = uLastIter;
        uCurrentIter -= deltaU;

        // do not clamp anything after 5 iterations
        static const int maxNumChoppedUpdates = getParam<int>("Newton.NumChoppedUpdates", 4);
        if (maxNumChoppedUpdates && this->numSteps_ <= maxNumChoppedUpdates)
        {
            lastIterationWasChopped_ = true;
            // clamp saturation change to at most 20% per iteration
            const auto& gridGeometry = this->assemblerReg().gridGeometry();
            auto fvGeometry = localView(gridGeometry);
            for (const auto& element : elements(gridGeometry.gridView()))
            {
                fvGeometry.bindElement(element);

                for (auto&& scv : scvs(fvGeometry))
                {
                    auto dofIdxGlobal = scv.dofIndex();

                    // calculate the old wetting phase saturation
                    const auto& spatialParams = this->assemblerReg().problem().spatialParams();
                    const auto elemSol = elementSolution(element, uCurrentIter, gridGeometry);

                    const auto fluidMatrixInteraction = spatialParams.fluidMatrixInteraction(element, scv, elemSol);
                    const double pw = uLastIter[dofIdxGlobal][0];
                    const double SnOld = uLastIter[dofIdxGlobal][1];

                    // const double SwOld = 1.0 - SnOld;
                    // const double pcOld = fluidMatrixInteraction.pc(SwOld);
                    // const double pn = pcOld + pw;
                    // const double pwMin = pn - fluidMatrixInteraction.pc(SwOld - 0.1);
                    // const double pwMax = pn - fluidMatrixInteraction.pc(SwOld + 0.1);

                    const double pwMin = pw - 1000;
                    const double pwMax = pw + 1000;

                    // clamp the result
                    using std::clamp;
                    uCurrentIter[dofIdxGlobal][0] = clamp(uCurrentIter[dofIdxGlobal][0], pwMin, pwMax);
                    uCurrentIter[dofIdxGlobal][1] = clamp(uCurrentIter[dofIdxGlobal][1], (SnOld - 0.1), (SnOld + 0.1));
                }
            }
        }
        else
            lastIterationWasChopped_ = false;

        // update the variables
        this->solutionChanged_(varsCurrentIter, uCurrentIter);

        if (this->enableResidualCriterion())
            this->computeResidualReduction_(varsCurrentIter);
    }

    bool lastIterationWasChopped_ = false;
    bool switchedInLastIteration_{false};
    std::vector<bool> preIndex_; // return the pore indicies which are predicted to be invaded
};

} // end namespace Dumux::PoreNetwork

#endif
