// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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
 * \ingroup PorousmediumCompositional
 * \copydoc Dumux::PriVarSwitchNewtonSolver
 */
#ifndef DUMUX_PRIVARSWITCH_NEWTON_SOLVER_HH
#define DUMUX_PRIVARSWITCH_NEWTON_SOLVER_HH

#include <memory>

#include <dumux/common/parameters.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumCompositional
 * \brief A newton solver that handles primary variable switches
 */
template <class Assembler, class LinearSolver, class PrimaryVariableSwitch>
class PriVarSwitchNewtonSolver : public NewtonSolver<Assembler, LinearSolver>
{
    using Scalar =  typename Assembler::Scalar;
    using ParentType = NewtonSolver<Assembler, LinearSolver>;
    using SolutionVector = typename Assembler::ResidualType;

public:
    using ParentType::ParentType;

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool newtonConverged() const override
    {
        if (switchedInLastIteration_)
            return false;

        return ParentType::newtonConverged();
    }

    /*!
     *
     * \brief Called before the Newton method is applied to an
     *        non-linear system of equations.
     *
     * \param u The initial solution
     */
    void newtonBegin(const SolutionVector &u) override
    {
        ParentType::newtonBegin(u);
        const int verbosity = getParamFromGroup<int>(this->paramGroup(), "PrimaryVariableSwitch.Verbosity", 1);
        priVarSwitch_ = std::make_unique<PrimaryVariableSwitch>(u.size(), verbosity);
    }

    /*!
     * \brief Indicates that one Newton iteration was finished.
     *
     * \param assembler The jacobian assembler
     * \param uCurrentIter The solution after the current Newton iteration
     * \param uLastIter The solution at the beginning of the current Newton iteration
     */
    void newtonEndStep(SolutionVector &uCurrentIter,
                       const SolutionVector &uLastIter) override
    {
        ParentType::newtonEndStep(uCurrentIter, uLastIter);

        // update the variable switch (returns true if the pri vars at at least one dof were switched)
        // for disabled grid variable caching
        auto& assembler = this->assembler();
        const auto& fvGridGeometry = assembler.fvGridGeometry();
        const auto& problem = assembler.problem();
        auto& gridVariables = assembler.gridVariables();

        // invoke the primary variable switch
        switchedInLastIteration_ = priVarSwitch_->update(uCurrentIter, gridVariables,
                                                         problem, fvGridGeometry);

        if(switchedInLastIteration_)
        {
            for (const auto& element : elements(fvGridGeometry.gridView()))
            {
                // if the volume variables are cached globally, we need to update those where the primary variables have been switched
                priVarSwitch_->updateSwitchedVolVars(problem, element, fvGridGeometry, gridVariables, uCurrentIter);

                // if the flux variables are cached globally, we need to update those where the primary variables have been switched
                priVarSwitch_->updateSwitchedFluxVarsCache(problem, element, fvGridGeometry, gridVariables, uCurrentIter);
            }
        }
    }

    /*!
     * \brief Called if the Newton method ended
     *        (not known yet if we failed or succeeded)
     */
    void newtonEnd() override
    {
        ParentType::newtonEnd();

        // in any way, we have to reset the switch flag
        switchedInLastIteration_ = false;
        // free some memory
        priVarSwitch_.release();
    }

private:
    //! the class handling the primary variable switch
    std::unique_ptr<PrimaryVariableSwitch> priVarSwitch_;
    //! if we switched primary variables in the last iteration
    bool switchedInLastIteration_ = false;
};

} // end namespace Dumux

#endif
