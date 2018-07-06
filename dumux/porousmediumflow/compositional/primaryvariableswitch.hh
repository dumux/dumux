// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
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
 * \brief The primary variable switch base class for compositional models
 */
#ifndef DUMUX_PRIMARY_VARIABLE_SWITCH_HH
#define DUMUX_PRIMARY_VARIABLE_SWITCH_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux {

/*!
 * \ingroup ImplicitModel
 * \brief Empty class for models without pri var switch
 */
class NoPrimaryVariableSwitch
{
public:
    NoPrimaryVariableSwitch(...) {}
    void init(...) {}
    bool wasSwitched(...) const { return false; }
    bool update(...) { return false; }
    void updateSwitchedVolVars(...) {}
    void updateSwitchedFluxVarsCache(...) {}
    bool update_(...) {return false; }
};

/*!
 * \ingroup PorousmediumCompositional
 * \brief The primary variable switch controlling the phase presence state variable
 */
template<class Implementation>
class PrimaryVariableSwitch
{
public:
    PrimaryVariableSwitch(const std::size_t& numDofs)
    {
        wasSwitched_.resize(numDofs, false);
    }

    //! If the primary variables were recently switched
    bool wasSwitched(std::size_t dofIdxGlobal) const
    {
        return wasSwitched_[dofIdxGlobal];
    }

    /*!
     * \brief Update the variable switch / phase presence
     *
     * \param curSol The current solution to be updated / modified
     * \param gridVariables The secondary variables on the grid
     * \param problem The problem
     * \param fvGridGeometry The finite-volume grid geometry
     */
    template<class SolutionVector, class GridVariables, class Problem>
    bool update(SolutionVector& curSol,
                GridVariables& gridVariables,
                const Problem& problem,
                const typename GridVariables::GridGeometry& fvGridGeometry)
    {
        bool switched = false;
        visited_.assign(wasSwitched_.size(), false);
        for (const auto& element : elements(fvGridGeometry.gridView()))
        {
            // make sure FVElementGeometry is bound to the element
            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            const auto curElemSol = elementSolution(element, curSol, fvGridGeometry);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto dofIdxGlobal = scv.dofIndex();
                if (!visited_[dofIdxGlobal])
                {
                    // Note this implies that volume variables don't differ
                    // in any sub control volume associated with the dof!
                    visited_[dofIdxGlobal] = true;
                    // Compute volVars on which grounds we decide
                    // if we need to switch the primary variables
                    auto& volVars = getVolVarAccess(gridVariables.curGridVolVars(), elemVolVars, scv);
                    volVars.update(curElemSol, problem, element, scv);

                    if (asImp_().update_(curSol[dofIdxGlobal], volVars, dofIdxGlobal, scv.dofPosition()))
                        switched = true;
                }
            }
        }

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag for our partition.
        if (fvGridGeometry.gridView().comm().size() > 1)
            switched = fvGridGeometry.gridView().comm().max(switched);

        return switched;
    }

    /*!
     * \brief Update the volume variables whose primary variables were
              switched. Required when volume variables are cached globally.
     */
    template<class Problem, class GridVariables, class SolutionVector,
             std::enable_if_t<GridVariables::GridVolumeVariables::cachingEnabled, int> = 0>
    void updateSwitchedVolVars(const Problem& problem,
                               const typename GridVariables::GridGeometry::GridView::template Codim<0>::Entity& element,
                               const typename GridVariables::GridGeometry& fvGridGeometry,
                               GridVariables& gridVariables,
                               const SolutionVector& sol)
    {
        // make sure FVElementGeometry is bound to the element
        auto fvGeometry = localView(fvGridGeometry);
        fvGeometry.bindElement(element);

        // update the secondary variables if global caching is enabled
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto dofIdxGlobal = scv.dofIndex();
            if (asImp_().wasSwitched(dofIdxGlobal))
            {
                const auto elemSol = elementSolution(element, sol, fvGridGeometry);
                auto& volVars = gridVariables.curGridVolVars().volVars(scv);
                volVars.update(elemSol, problem, element, scv);
            }
        }
    }

    /*!
     * \brief Update the fluxVars cache for dof whose primary variables were
              switched. Required when flux variables are cached globally (not for box method).
     */
     template<class Problem, class GridVariables, class SolutionVector,
              std::enable_if_t<(GridVariables::GridFluxVariablesCache::cachingEnabled &&
                                GridVariables::GridGeometry::discMethod != DiscretizationMethod::box), int> = 0>
     void updateSwitchedFluxVarsCache(const Problem& problem,
                                const typename GridVariables::GridGeometry::GridView::template Codim<0>::Entity& element,
                                const typename GridVariables::GridGeometry& fvGridGeometry,
                                GridVariables& gridVariables,
                                const SolutionVector& sol)
    {
        // update the flux variables if global caching is enabled
        const auto dofIdxGlobal = fvGridGeometry.dofMapper().index(element);

        if (asImp_().wasSwitched(dofIdxGlobal))
        {
            // make sure FVElementGeometry and the volume variables are bound
            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bind(element);
            auto curElemVolVars = localView(gridVariables.curGridVolVars());
            curElemVolVars.bind(element, fvGeometry, sol);
            gridVariables.gridFluxVarsCache().updateElement(element, fvGeometry, curElemVolVars);
        }
    }

     //! brief Do nothing when volume variables are not cached globally.
     template<class Problem, class GridVariables, class SolutionVector,
              std::enable_if_t<!GridVariables::GridVolumeVariables::cachingEnabled, int> = 0>
    void updateSwitchedVolVars(const Problem& problem,
                               const typename GridVariables::GridGeometry::GridView::template Codim<0>::Entity& element,
                               const typename GridVariables::GridGeometry& fvGridGeometry,
                               GridVariables& gridVariables,
                               const SolutionVector &uCurrentIter) const {}

    //! brief Do nothing when flux variables are not cached globally or the box method is used.
    template<class Problem, class GridVariables, class SolutionVector,
             std::enable_if_t<(!GridVariables::GridFluxVariablesCache::cachingEnabled ||
                               GridVariables::GridGeometry::discMethod == DiscretizationMethod::box), int> = 0>
    void updateSwitchedFluxVarsCache(const Problem& problem,
                               const typename GridVariables::GridGeometry::GridView::template Codim<0>::Entity& element,
                               const typename GridVariables::GridGeometry& fvGridGeometry,
                               GridVariables& gridVariables,
                               const SolutionVector& sol) const {}

protected:

    //! return actual implementation (static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    //! return actual implementation (static polymorphism)
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    // perform variable switch at a degree of freedom location
    template<class VolumeVariables, class GlobalPosition>
    bool update_(typename VolumeVariables::PrimaryVariables& priVars,
                 const VolumeVariables& volVars,
                 std::size_t dofIdxGlobal,
                 const GlobalPosition& globalPos)
    {
        // evaluate if the primary variable switch would switch
        // to be implemented by the deriving class
        DUNE_THROW(Dune::NotImplemented, "This model seems to use a primary variable switch but none is implemented!");
    }

    std::vector<bool> wasSwitched_;
    std::vector<bool> visited_;

private:
    template<class GridVolumeVariables, class ElementVolumeVariables, class SubControlVolume>
    static auto getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    -> std::enable_if_t<!GridVolumeVariables::cachingEnabled, decltype(elemVolVars[scv])>
    { return elemVolVars[scv]; }

    template<class GridVolumeVariables, class ElementVolumeVariables, class SubControlVolume>
    static auto getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    -> std::enable_if_t<GridVolumeVariables::cachingEnabled, decltype(gridVolVars.volVars(scv))>
    { return gridVolVars.volVars(scv); }
};

} // end namespace dumux

#endif
