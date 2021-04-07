// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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
 * \ingroup Nonlinear
 * \ingroup MultiDomain
 * \copydoc Dumux::MultiDomainNewtonSolver
 */
#ifndef DUMUX_MULTIDOMAIN_NEWTON_SOLVER_HH
#define DUMUX_MULTIDOMAIN_NEWTON_SOLVER_HH

#include <memory>
#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux {
namespace Detail {

template<class Assembler, class Index>
using DetectPVSwitchMultiDomain = typename Assembler::template GridVariables<Index::value>::VolumeVariables::PrimaryVariableSwitch;

template<class Assembler, std::size_t i>
using GetPVSwitchMultiDomain = Dune::Std::detected_or<int, DetectPVSwitchMultiDomain, Assembler, Dune::index_constant<i>>;

} // end namespace Detail

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
    using typename ParentType::Backend;
    using typename ParentType::SolutionVector;

    static constexpr bool assemblerExportsVariables = Detail::exportsVariables<Assembler>;

    template<std::size_t i>
    using PrimaryVariableSwitch = typename Detail::GetPVSwitchMultiDomain<Assembler, i>::type;

    template<std::size_t i>
    using HasPriVarsSwitch = typename Detail::GetPVSwitchMultiDomain<Assembler, i>::value_t; // std::true_type or std::false_type

    template<std::size_t i>
    using PrivarSwitchPtr = std::unique_ptr<PrimaryVariableSwitch<i>>;
    using PriVarSwitchPtrTuple = typename Assembler::Traits::template Tuple<PrivarSwitchPtr>;

public:
    using typename ParentType::Variables;

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
    void newtonBeginStep(const Variables& varsCurrentIter) override
    {
        ParentType::newtonBeginStep(varsCurrentIter);
        couplingManager_->updateSolution(Backend::dofs(varsCurrentIter));
    }

    /*!
     *
     * \brief Called before the Newton method is applied to an
     *        non-linear system of equations.
     *
     * \param vars The variables representing the initial solution
     */
    void newtonBegin(Variables& vars) override
    {
        ParentType::newtonBegin(vars);

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Assembler::Traits::numSubDomains>{}, [&](auto&& id)
        {
            this->initPriVarSwitch_(vars, id, HasPriVarsSwitch<std::decay_t<decltype(id)>::value>{});
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
     * \param varsCurrentIter The variables after the current Newton iteration
     * \param uLastIter The solution at the beginning of the current Newton iteration
     */
    void newtonEndStep(Variables& varsCurrentIter, const SolutionVector& uLastIter) override
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Assembler::Traits::numSubDomains>{}, [&](auto&& id)
        {
            auto& uCurrentIter = Backend::dofs(varsCurrentIter)[id];
            if constexpr (!assemblerExportsVariables)
                this->invokePriVarSwitch_(this->assembler().gridVariables(id),
                                          uCurrentIter, id, HasPriVarsSwitch<std::decay_t<decltype(id)>::value>{});
            else
                this->invokePriVarSwitch_(varsCurrentIter[id], uCurrentIter, id, HasPriVarsSwitch<std::decay_t<decltype(id)>::value>{});
        });

        ParentType::newtonEndStep(varsCurrentIter, uLastIter);
        couplingManager_->updateSolution(Backend::dofs(varsCurrentIter));
    }

protected:
    /*!
     * \brief Update solution-depended quantities like grid variables after the solution has changed.
     */
    void solutionChanged_(Variables& vars, const SolutionVector& uCurrentIter) override
    {
        couplingManager_->updateSolution(uCurrentIter);
        ParentType::solutionChanged_(vars, uCurrentIter);
    }

private:

    /*!
     * \brief Reset the privar switch state, noop if there is no priVarSwitch
     */
    template<std::size_t i>
    void initPriVarSwitch_(Variables&, Dune::index_constant<i> id, std::false_type) {}

    /*!
     * \brief Switch primary variables if necessary
     */
    template<std::size_t i>
    void initPriVarSwitch_(Variables& vars, Dune::index_constant<i> id, std::true_type)
    {
        using namespace Dune::Hybrid;
        auto& priVarSwitch = *elementAt(priVarSwitches_, id);
        auto& sol = Backend::dofs(vars)[id];

        priVarSwitch.reset(sol.size());
        priVarsSwitchedInLastIteration_[i] = false;

        const auto& problem = this->assembler().problem(id);
        const auto& gridGeometry = this->assembler().gridGeometry(id);
        if constexpr (!assemblerExportsVariables)
            priVarSwitch.updateDirichletConstraints(problem, gridGeometry, this->assembler().gridVariables(id), sol);
        else // This expects usage of MultiDomainGridVariables
            priVarSwitch.updateDirichletConstraints(problem, gridGeometry, vars[id], sol[id]);
    }

    /*!
     * \brief Switch primary variables if necessary, noop if there is no priVarSwitch
     */
    template<class SubVars, class SubSol, std::size_t i>
    void invokePriVarSwitch_(SubVars&, SubSol&, Dune::index_constant<i> id, std::false_type) {}

    /*!
     * \brief Switch primary variables if necessary
     */
    template<class SubVars, class SubSol, std::size_t i>
    void invokePriVarSwitch_(SubVars& subVars, SubSol& uCurrentIter, Dune::index_constant<i> id, std::true_type)
    {
        // update the variable switch (returns true if the pri vars at at least one dof were switched)
        // for disabled grid variable caching
        const auto& gridGeometry = this->assembler().gridGeometry(id);
        const auto& problem = this->assembler().problem(id);

        using namespace Dune::Hybrid;
        auto& priVarSwitch = *elementAt(priVarSwitches_, id);

        // invoke the primary variable switch
        priVarsSwitchedInLastIteration_[i] = priVarSwitch.update(uCurrentIter, subVars, problem, gridGeometry);

        if (priVarsSwitchedInLastIteration_[i])
        {
            for (const auto& element : elements(gridGeometry.gridView()))
            {
                // if the volume variables are cached globally, we need to update those where the primary variables have been switched
                priVarSwitch.updateSwitchedVolVars(problem, element, gridGeometry, subVars, uCurrentIter);

                // if the flux variables are cached globally, we need to update those where the primary variables have been switched
                priVarSwitch.updateSwitchedFluxVarsCache(problem, element, gridGeometry, subVars, uCurrentIter);
            }
        }
    }

    //! the coupling manager
    std::shared_ptr<CouplingManager> couplingManager_;

    //! the class handling the primary variable switch
    PriVarSwitchPtrTuple priVarSwitches_;
    //! if we switched primary variables in the last iteration
    std::array<bool, Assembler::Traits::numSubDomains> priVarsSwitchedInLastIteration_;
};

} // end namespace Dumux

#endif
