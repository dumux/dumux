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
 * \ingroup TimeStepping
 * \brief A time stepper performing a single time step of a transient simulation
 */
#ifndef DUMUX_TIMESTEPPING_MULTISTAGE_TIMESTEPPER_HH
#define DUMUX_TIMESTEPPING_MULTISTAGE_TIMESTEPPER_HH

#include <memory>
#include <type_traits>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/common/pdesolver.hh>
#include <dumux/experimental/new_assembly/dumux/common/variables.hh>

#include <dumux/experimental/new_assembly/dumux/timestepping/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/timestepping/multistageparams.hh>
#include <dumux/experimental/new_assembly/dumux/timestepping/multistagemethods.hh>

namespace Dumux {

/*!
 * \ingroup TimeStepping
 * \brief Time stepping with a multi-stage method
 * \tparam V The type used for variables of the PDE system
 * \note We limit ourselves to "diagonally" implicit multi-stage methods where solving
 *       a stage can only depend on the values of the same stage and stages before
 *       but not future stages (which would require solving larger linear systems)
 */
template<Concepts::MultiStageVariables V>
class MultiStageTimeStepper
{
    using Scalar = Variables::ScalarType<V>;
    using StageParams = MultiStageParams<Scalar>;

public:
    using Variables = V;
    using StageVariables = typename V::StageVariables;

    /*!
     * \brief The constructor
     * \param pdeSolver Solver class for solving a PDE in each stage
     * \param msMethod The multi-stage method which is to be used for time integration
     * \todo TODO: Add time step control if the pde solver doesn't converge
     */
    MultiStageTimeStepper(std::shared_ptr<PDESolver<Variables>> pdeSolver,
                          std::shared_ptr<const MultiStageMethod<Scalar>> msMethod)
    : pdeSolver_(pdeSolver)
    , msMethod_(msMethod)
    {}

    /*!
     * \brief Advance one time step of the given time loop
     * \param vars The variables object at the current time level.
     * \param t The current time level
     * \param dt The time step size to be used
     * \note We expect the time level in vars to correspond to the given time `t`
     * \todo: TODO: Add time step control if the pde solver doesn't converge
     */
    void step(Variables& vars, const Scalar t, const Scalar dt)
    {
        vars.clearStages();
        for (auto stageIdx = 1UL; stageIdx <= msMethod_->numStages(); ++stageIdx)
            solveStage_(vars, std::make_shared<StageParams>(*msMethod_, stageIdx, t, dt));
        vars.clearStages();
    }

    /*!
     * \brief Overload for given single-stage variables, which wraps
     *        the stage variables into `MultiStageVariables`. This
     *        overload is only visible if the given variables are a
     *        view on stage variables that can be constructed from a reference.
     */
    void step(StageVariables& vars, const Scalar t, const Scalar dt) requires(
        Concepts::View<V> and std::constructible_from<Variables, StageVariables&>)
    {
        Variables multiStageVars{vars};
        step(multiStageVars, t, dt);
    }

    //! Set/change the time step method
    void setMethod(std::shared_ptr<const MultiStageMethod<Scalar>> msMethod)
    { msMethod_ = msMethod; }

private:
    void solveStage_(Variables& vars, std::shared_ptr<StageParams> params)
    {
        vars.setStageParams(params);
        pdeSolver_->solve(vars);
    }

    std::shared_ptr<PDESolver<Variables>> pdeSolver_;
    std::shared_ptr<const MultiStageMethod<Scalar>> msMethod_;
};

template<typename Solver, typename Method> requires(
    std::derived_from<Solver, PDESolver<typename Solver::Variables>>)
MultiStageTimeStepper(std::shared_ptr<Solver>,
                      std::shared_ptr<Method>) -> MultiStageTimeStepper<typename Solver::Variables>;

} // namespace Dumux

#endif
