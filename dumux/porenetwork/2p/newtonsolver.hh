// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMTwoPModel
 * \copydoc Dumux::PoreNetwork::TwoPNewtonSolver
 */
#ifndef DUMUX_PNM_NEWTON_SOLVER_HH
#define DUMUX_PNM_NEWTON_SOLVER_HH

#include <memory>
#include <string>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/porenetwork/2p/invasionstate.hh>
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
//! Forward declaration of the two-phase PNM newton solver, specialized on the state switch method
template<class Assembler, class LinearSolver,
         template<class, class> class NewtonConsistencyChecks = TwoPNewtonConsistencyChecks,
         class SwitchMethod = StateSwitchMethod::IterationWithPrevTimeStepInformation>
class TwoPNewtonSolver;

template<class Assembler, class LinearSolver,
         template<class, class> class NewtonConsistencyChecks>
class TwoPNewtonSolver<Assembler, LinearSolver, NewtonConsistencyChecks,
                       StateSwitchMethod::IterationWithPrevTimeStepInformation>
: public Dumux::NewtonSolver<Assembler, LinearSolver>
{
    using ParentType =  Dumux::NewtonSolver<Assembler, LinearSolver>;
    using SolutionVector = typename ParentType::SolutionVector;

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
     *        be accurate enough. We enforce an additional Newton iteration if the invasion
     *        switch was triggered in the last iteration.
     */
    bool newtonConverged() const final
    {
        if (switchedInLastIteration_)
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

private:
    bool switchedInLastIteration_{false};
};

/*!
 * \ingroup PNMTwoPModel
 * \brief A two-phase PNM specific newton solver switching the invasion state per iteration.
 *
 * In contrast to the deprecated form, the invasion states of the current and the previous time
 * level are handed to the solver, which is the only place that changes them during the time step.
 */
template<class Assembler, class LinearSolver,
         template<class, class> class NewtonConsistencyChecks>
class TwoPNewtonSolver<Assembler, LinearSolver, NewtonConsistencyChecks, StateSwitchMethod::Iteration>
: public Dumux::NewtonSolver<Assembler, LinearSolver>
{
    using ParentType = Dumux::NewtonSolver<Assembler, LinearSolver>;
    using SolutionVector = typename ParentType::SolutionVector;
    using InvasionState = typename Assembler::Problem::InvasionState;

public:
    /*!
     * \brief The constructor
     * \param assembler The assembler
     * \param linearSolver The linear solver
     * \param invasionState The invasion state of the current time level
     * \param prevInvasionState The invasion state of the previous time level
     * \param comm The communicator
     * \param paramGroup The parameter group
     */
    TwoPNewtonSolver(std::shared_ptr<Assembler> assembler,
                     std::shared_ptr<LinearSolver> linearSolver,
                     std::shared_ptr<InvasionState> invasionState,
                     std::shared_ptr<InvasionState> prevInvasionState,
                     const typename ParentType::Communication& comm = Dune::MPIHelper::getCommunication(),
                     const std::string& paramGroup = "")
    : ParentType(assembler, linearSolver, comm, paramGroup)
    , invasionState_(invasionState)
    , prevInvasionState_(prevInvasionState)
    {}

    //! Called after each Newton update
    void newtonEndStep(SolutionVector& uCurrentIter, const SolutionVector& uLastIter) final
    {
        ParentType::newtonEndStep(uCurrentIter, uLastIter);

        auto& gridVariables = this->assembler().gridVariables();
        switchedInLastIteration_ = invasionState_->update(
            uCurrentIter, gridVariables.curGridVolVars(), gridVariables.gridFluxVarsCache()
        );

        // If the solution is about to be accepted, check for accuracy and trigger a retry
        // with a decreased time step size if necessary.
        if (newtonConverged())
        {
            NewtonConsistencyChecks<typename Assembler::GridVariables, SolutionVector> checks;
            checks.performChecks(gridVariables, uCurrentIter, this->assembler().prevSol(),
                                 *invasionState_, *prevInvasionState_);
        }
    }

    //! We enforce an additional Newton iteration if the invasion switch was triggered
    bool newtonConverged() const final
    {
        if (switchedInLastIteration_)
            return false;

        return ParentType::newtonConverged();
    }

    //! Restore the invasion state of the previous time level
    void newtonFail(SolutionVector& u) final
    {
        ParentType::newtonFail(u);
        *invasionState_ = *prevInvasionState_;
    }

private:
    std::shared_ptr<InvasionState> invasionState_;
    std::shared_ptr<InvasionState> prevInvasionState_;
    bool switchedInLastIteration_{false};
};

} // end namespace Dumux::PoreNetwork

#endif
