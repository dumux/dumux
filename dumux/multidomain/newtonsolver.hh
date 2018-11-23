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
 * \copydoc Dumux::MultiDomainNewtonSolver
 */
#ifndef DUMUX_MULTIDOMAIN_NEWTON_SOLVER_HH
#define DUMUX_MULTIDOMAIN_NEWTON_SOLVER_HH

#include <memory>
#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux {

namespace MDDetail {

template<class Assembler, class IdType>
using GridVariables = std::decay_t<decltype(std::declval<Assembler>().gridVariables(IdType()))>;

template<class Assembler, class IdType>
using GetPVSwitch = typename GridVariables<Assembler, IdType>::VolumeVariables::PrimaryVariableSwitch;

template<class Assembler, std::size_t i>
using PrimaryVariableSwitch = Dune::Std::detected_or<int, GetPVSwitch, Assembler, Dune::index_constant<i>>;

}

/*!
 * \ingroup Nonlinear
 * \ingroup MultiDomain
 * \brief Newton sover for coupled problems
 */
template <class Assembler, class LinearSolver, class CouplingManager,
          class Reassembler = DefaultPartialReassembler,
          class Comm = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> >
class MultiDomainNewtonSolver: public NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>
{
    using ParentType = NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>;
    using SolutionVector = typename Assembler::ResidualType;

    template<std::size_t i>
    using PrimaryVariableSwitch = typename MDDetail::PrimaryVariableSwitch<Assembler, i>::type;

    template<std::size_t i>
    using PrivarSwitchPtr = std::unique_ptr<PrimaryVariableSwitch<i>>;

    using PriVarSwitchPtrTuple = typename Assembler::Traits::template MakeTuple<PrivarSwitchPtr>;

    template<std::size_t i>
    using HasPriVarsSwitch = typename MDDetail::PrimaryVariableSwitch<Assembler, i>::value_t;

public:

    /*!
     * \brief The constructor
     */
    MultiDomainNewtonSolver(std::shared_ptr<Assembler> assembler,
                            std::shared_ptr<LinearSolver> linearSolver,
                            std::shared_ptr<CouplingManager> couplingManager,
                            const Comm& comm = Dune::MPIHelper::getCollectiveCommunication(),
                            const std::string& paramGroup = "")
    : ParentType(assembler, linearSolver, comm, paramGroup)
    , couplingManager_(couplingManager)
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Assembler::Traits::numSubDomains>{}, [&](auto&& id)
        {
            const int priVarSwitchVerbosity = getParamFromGroup<int>(paramGroup, "PrimaryVariableSwitch.Verbosity", 1);
            using PVSwitch = PrimaryVariableSwitch<std::decay_t<decltype(id)>::value>;
            elementAt(priVarSwitches_, id) = std::make_unique<PVSwitch>(priVarSwitchVerbosity);
        });

        priVarsSwitchedInLastIteration_.fill(false);
    }

    /*!
     * \brief Indicates the beginning of a Newton iteration.
     */
    void newtonBeginStep(const SolutionVector& uCurrentIter) override
    {
        ParentType::newtonBeginStep(uCurrentIter);
        couplingManager_->updateSolution(uCurrentIter);
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

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Assembler::Traits::numSubDomains>{}, [&](auto&& id)
        {
            resetPriVarSwitch_(u[id].size(), id, HasPriVarsSwitch<id>());
        });
    }

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool newtonConverged() const override
    {
        if (Dune::any_true(priVarsSwitchedInLastIteration_))
            return false;

        return ParentType::newtonConverged();
    }


    /*!
     * \brief Indicates that one Newton iteration was finished.
     *
     * \param assembler The jacobian assembler
     * \param uCurrentIter The solution after the current Newton iteration
     * \param uLastIter The solution at the beginning of the current Newton iteration
     */
    void newtonEndStep(SolutionVector &uCurrentIter, const SolutionVector &uLastIter) override
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Assembler::Traits::numSubDomains>{}, [&](auto&& id)
        {
            invokePriVarSwitch_(uCurrentIter[id], id, HasPriVarsSwitch<id>());
        });

        ParentType::newtonEndStep(uCurrentIter, uLastIter);
        couplingManager_->updateSolution(uCurrentIter);
    }

private:


    template<std::size_t i>
    void resetPriVarSwitch_(const std::size_t numDofs, Dune::index_constant<i> id, std::false_type)
    {}

    template<std::size_t i>
    void resetPriVarSwitch_(const std::size_t numDofs, Dune::index_constant<i> id, std::true_type)
    {
        using namespace Dune::Hybrid;
        elementAt(priVarSwitches_, id)->reset(numDofs);
        priVarsSwitchedInLastIteration_[i] = false;
    }

    template<class SubSol, std::size_t i>
    void invokePriVarSwitch_(SubSol&, Dune::index_constant<i> id, std::false_type)
    {}

    template<class SubSol, std::size_t i>
    void invokePriVarSwitch_(SubSol& uCurrentIter, Dune::index_constant<i> id, std::true_type)
    {
        // update the variable switch (returns true if the pri vars at at least one dof were switched)
        // for disabled grid variable caching
        const auto& fvGridGeometry = this->assembler().fvGridGeometry(id);
        const auto& problem = this->assembler().problem(id);
        auto& gridVariables = this->assembler().gridVariables(id);

        using namespace Dune::Hybrid;
        auto& priVarSwitch = *elementAt(priVarSwitches_, id);

        // invoke the primary variable switch
        priVarsSwitchedInLastIteration_[i] = priVarSwitch.update(uCurrentIter, gridVariables,
                                                                 problem, fvGridGeometry);

        if (priVarsSwitchedInLastIteration_[i])
        {
            for (const auto& element : elements(fvGridGeometry.gridView()))
            {
                // if the volume variables are cached globally, we need to update those where the primary variables have been switched
                priVarSwitch.updateSwitchedVolVars(problem, element, fvGridGeometry, gridVariables, uCurrentIter);

                // if the flux variables are cached globally, we need to update those where the primary variables have been switched
                priVarSwitch.updateSwitchedFluxVarsCache(problem, element, fvGridGeometry, gridVariables, uCurrentIter);
            }
        }
    }


    std::shared_ptr<CouplingManager> couplingManager_;

    //! the class handling the primary variable switch
    PriVarSwitchPtrTuple priVarSwitches_;
    //! if we switched primary variables in the last iteration
    std::array<bool, Assembler::Traits::numSubDomains> priVarsSwitchedInLastIteration_;
};

} // end namespace Dumux

#endif
