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
 * \brief An adapter for the Newton to manage models with primary variable switch
 */
#ifndef DUMUX_NONLINEAR_PRIMARY_VARIABLE_SWITCH_ADAPTER_HH
#define DUMUX_NONLINEAR_PRIMARY_VARIABLE_SWITCH_ADAPTER_HH

#include <memory>
#include <dune/common/std/type_traits.hh>

namespace Dumux {
namespace Impl {

//! helper aliases to extract a primary variable switch from the VolumeVariables (if defined, yields int otherwise)
template<class Variables>
using DetectPVSwitch = typename Variables::VolumeVariables::PrimaryVariableSwitch;

template<class Variables>
using GetPVSwitch = Dune::Std::detected_or<int, DetectPVSwitch, Variables>;

template<class Variables>
using PrimaryVariableSwitch = typename GetPVSwitch<Variables>::type;

template<class Variables>
constexpr bool hasPriVarsSwitch() { return typename GetPVSwitch<Variables>::value_t(); };

} // end namespace Impl

/*!
 * \ingroup Nonlinear
 * \brief An adapter for the Newton to manage models with primary variable switch
 */
template <class Variables, bool isValid = Impl::hasPriVarsSwitch<Variables>()>
class PrimaryVariableSwitchAdapter
{
    using PrimaryVariableSwitch = typename Impl::PrimaryVariableSwitch<Variables>;

public:
    PrimaryVariableSwitchAdapter(const std::string& paramGroup = "")
    {
        const int priVarSwitchVerbosity = getParamFromGroup<int>(paramGroup, "PrimaryVariableSwitch.Verbosity", 1);
        priVarSwitch_ = std::make_unique<PrimaryVariableSwitch>(priVarSwitchVerbosity);
    }

    /*!
     * \brief Initialize the privar switch
     */
    template<class SolutionVector>
    void initialize(SolutionVector& sol, Variables& vars)
    {
        priVarSwitch_->reset(sol.size());
        priVarsSwitchedInLastIteration_ = false;
        const auto& problem = vars.curGridVolVars().problem();
        const auto& gridGeometry = problem.gridGeometry();
        priVarSwitch_->updateBoundary(problem, gridGeometry, vars, sol);
    }

    /*!
     * \brief Switch primary variables if necessary
     */
    template<class SolutionVector>
    void invoke(SolutionVector& uCurrentIter, Variables& vars)
    {
        // update the variable switch (returns true if the pri vars at at least one dof were switched)
        // for disabled grid variable caching
        const auto& problem = vars.curGridVolVars().problem();
        const auto& gridGeometry = problem.gridGeometry();

        // invoke the primary variable switch
        priVarsSwitchedInLastIteration_ = priVarSwitch_->update(uCurrentIter, vars, problem, gridGeometry);
        if (priVarsSwitchedInLastIteration_)
        {
            for (const auto& element : elements(gridGeometry.gridView()))
            {
                // if the volume variables are cached globally, we need to update those where the primary variables have been switched
                priVarSwitch_->updateSwitchedVolVars(problem, element, gridGeometry, vars, uCurrentIter);

                // if the flux variables are cached globally, we need to update those where the primary variables have been switched
                priVarSwitch_->updateSwitchedFluxVarsCache(problem, element, gridGeometry, vars, uCurrentIter);
            }
        }
    }

    /*!
     * \brief Whether the primary variables have been switched in the last call to invoke
     */
    bool switched() const
    { return priVarsSwitchedInLastIteration_; }

private:
    //! the class handling the primary variable switch
    std::unique_ptr<PrimaryVariableSwitch> priVarSwitch_;
    //! if we switched primary variables in the last iteration
    bool priVarsSwitchedInLastIteration_ = false;
};

/*!
 * \ingroup Nonlinear
 * \brief An empty adapter for the Newton for models without primary variable switch
 */
template <class Variables>
class PrimaryVariableSwitchAdapter<Variables, false>
{
public:
    PrimaryVariableSwitchAdapter(const std::string& paramGroup = "") {}

    template<class SolutionVector>
    void initialize(SolutionVector&, Variables&) {}

    template<class SolutionVector>
    void invoke(SolutionVector&, Variables&) {}

    bool switched() const { return false; }
};

} // end namespace Dumux

#endif
