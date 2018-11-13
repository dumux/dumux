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
 * \ingroup Nonlinear
 * \ingroup MultiDomain
 * \copydoc Dumux::MultiDomainPriVarSwitchNewtonSolver
 */
#ifndef DUMUX_MULTIDOMAIN_PRIVARSWITCH_NEWTON_SOLVER_HH
#define DUMUX_MULTIDOMAIN_PRIVARSWITCH_NEWTON_SOLVER_HH

#include <memory>

#include <dumux/common/typetraits/utility.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/multidomain/newtonsolver.hh>

namespace Dumux {

/*!
 * \ingroup Nonlinear
 * \ingroup MultiDomain
 * \brief A newton solver that handles primary variable switches for multi domain problems
 */
template <class Assembler, class LinearSolver, class CouplingManager, class PrivarSwitchTypeTuple,
          class Reassembler = DefaultPartialReassembler,
          class Comm = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> >
class MultiDomainPriVarSwitchNewtonSolver : public MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager, Reassembler, Comm>
{
    using Scalar =  typename Assembler::Scalar;
    using ParentType = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager, Reassembler, Comm>;
    using SolutionVector = typename Assembler::ResidualType;

    template<std::size_t i> using PriVarSwitch = std::tuple_element_t<i, PrivarSwitchTypeTuple>;
    template<std::size_t i> using PrivarSwitchPtr = std::unique_ptr<PriVarSwitch<i>>;
    using PriVarSwitchPtrTuple = typename Assembler::Traits::template MakeTuple<PrivarSwitchPtr>;

public:
    using ParentType::ParentType;

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool newtonConverged() const override
    {
        if (std::any_of(switchedInLastIteration_.begin(), switchedInLastIteration_.end(), [](bool i){ return i;}))
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
        switchedInLastIteration_.fill(false);

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Assembler::Traits::numSubDomains>{}, [&](auto&& id)
        {
            using PVSwitch = PriVarSwitch<std::decay_t<decltype(id)>::value>;
            elementAt(priVarSwitches_, id) = std::make_unique<PVSwitch>(u[id].size());
        });
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

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Assembler::Traits::numSubDomains>{}, [&](auto&& id)
        {
            const auto& fvGridGeometry = assembler.fvGridGeometry(id);
            const auto& problem = assembler.problem(id);
            auto& gridVariables = assembler.gridVariables(id);
            auto& priVarSwitch = elementAt(priVarSwitches_, id);

            // invoke the primary variable switch
            switchedInLastIteration_[id] = priVarSwitch->update(uCurrentIter[id], gridVariables,
                                                                               problem, fvGridGeometry);

            if(switchedInLastIteration_[id])
            {
                for (const auto& element : elements(fvGridGeometry.gridView()))
                {
                    // if the volume variables are cached globally, we need to update those where the primary variables have been switched
                    priVarSwitch->updateSwitchedVolVars(problem, element, fvGridGeometry, gridVariables, uCurrentIter);

                    // if the flux variables are cached globally, we need to update those where the primary variables have been switched
                    priVarSwitch->updateSwitchedFluxVarsCache(problem, element, fvGridGeometry, gridVariables, uCurrentIter);
                }
            }
        });
    }

    /*!
     * \brief Called if the Newton method ended
     *        (not known yet if we failed or succeeded)
     */
    void newtonEnd() override
    {
        ParentType::newtonEnd();

        // in any way, we have to reset the switch flag
        switchedInLastIteration_.fill(false);

        // free some memory
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Assembler::Traits::numSubDomains>{}, [&](auto&& id)
        {
            elementAt(priVarSwitches_, id).release();
        });
    }

private:
    //! the class handling the primary variable switch
    PriVarSwitchPtrTuple priVarSwitches_;
    //! if we switched primary variables in the last iteration
    std::array<bool, Assembler::Traits::numSubDomains> switchedInLastIteration_;
};

} // end namespace Dumux

#endif
