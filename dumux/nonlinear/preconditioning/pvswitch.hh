// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Nonlinear
 * \brief Applying a model's primary variable switch to one subdomain at a time
 */
#ifndef DUMUX_NONLINEAR_PRECONDITIONING_PVSWITCH_HH
#define DUMUX_NONLINEAR_PRECONDITIONING_PVSWITCH_HH

#include <cstddef>
#include <memory>
#include <string>

#include <dune/grid/common/rangegenerators.hh>

#include <dumux/common/parameters.hh>
#include <dumux/nonlinear/primaryvariableswitchadapter.hh>

namespace Dumux::NonlinearPreconditioning {

/*!
 * \ingroup Nonlinear
 * \brief Owns a model's primary variable switch and applies it globally or to a subdomain
 *
 * For models without a switch every member compiles to nothing, so solvers can call into this
 * unconditionally.
 *
 * The switch decision at a degree of freedom depends only on that degree of freedom's own primary and
 * volume variables, so restricting the loop to a subdomain's elements restricts the effect exactly, with
 * no coordination between subdomains. What the caller must guarantee is ownership: a degree of freedom
 * may be switched by one subdomain only, otherwise two local solves would each rewrite the same state
 * from different frozen data.
 */
template<class GridVariables>
class SubdomainPrimaryVariableSwitch
{
    static constexpr bool hasSwitch = hasPriVarsSwitch<GridVariables>;
    using Switch = typename Detail::PrimaryVariableSwitch<GridVariables>;

public:
    explicit SubdomainPrimaryVariableSwitch(const std::string& paramGroup = "")
    {
        if constexpr (hasSwitch)
        {
            const int verbosity = getParamFromGroup<int>(paramGroup, "PrimaryVariableSwitch.Verbosity", 0);
            switch_ = std::make_unique<Switch>(verbosity);
        }
    }

    static constexpr bool active() { return hasSwitch; }

    //! establish the initial states and any Dirichlet constraints, over the whole grid
    template<class SolutionVector, class Problem, class GridGeometry>
    void initialize(SolutionVector& sol, GridVariables& gridVariables,
                    const Problem& problem, const GridGeometry& gridGeometry)
    {
        if constexpr (hasSwitch)
        {
            switch_->reset(sol.size());
            switch_->updateDirichletConstraints(problem, gridGeometry, gridVariables, sol);
        }
    }

    //! apply over the whole grid, as an ordinary Newton step would
    template<class SolutionVector, class Problem, class GridGeometry>
    bool applyGlobal(SolutionVector& sol, GridVariables& gridVariables,
                     const Problem& problem, const GridGeometry& gridGeometry)
    {
        if constexpr (hasSwitch)
        {
            const bool switched = switch_->update(sol, gridVariables, problem, gridGeometry);
            if (switched)
                for (const auto& element : elements(gridGeometry.gridView()))
                    refreshCaches_(element, sol, gridVariables, problem, gridGeometry);
            return switched;
        }
        else
            return false;
    }

    /*!
     * \brief Apply to the given elements only
     * \param elementIndices the elements owned by the subdomain; for cell-centred schemes these are its
     *        degrees of freedom, so nothing outside the subdomain can be touched
     */
    template<class SolutionVector, class Problem, class GridGeometry, class ElementIndices>
    bool applyToSubdomain(SolutionVector& sol,
                          GridVariables& gridVariables,
                          const Problem& problem,
                          const GridGeometry& gridGeometry,
                          const ElementIndices& elementIndices)
    {
        if constexpr (hasSwitch)
        {
            const bool switched = switch_->update(sol, gridVariables, problem, gridGeometry, elementIndices);
            if (switched)
                for (const auto eIdx : elementIndices)
                    refreshCaches_(gridGeometry.element(eIdx), sol, gridVariables, problem, gridGeometry);
            return switched;
        }
        else
            return false;
    }

    std::size_t numSwitched() const
    {
        if constexpr (hasSwitch)
            return switch_->numSwitched();
        else
            return 0;
    }

private:
    template<class Element, class SolutionVector, class Problem, class GridGeometry>
    void refreshCaches_(const Element& element, const SolutionVector& sol, GridVariables& gridVariables,
                        const Problem& problem, const GridGeometry& gridGeometry)
    {
        if constexpr (hasSwitch)
        {
            switch_->updateSwitchedVolVars(problem, element, gridGeometry, gridVariables, sol);
            switch_->updateSwitchedFluxVarsCache(problem, element, gridGeometry, gridVariables, sol);
        }
    }

    std::unique_ptr<Switch> switch_;
};

} // end namespace Dumux::NonlinearPreconditioning

#endif
