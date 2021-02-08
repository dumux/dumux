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
 * \ingroup Discretization
 * \brief The grid variable class for finite volume schemes,
 *        storing variables on scv and scvf (volume and flux variables)
 */
#ifndef DUMUX_FV_GRID_VARIABLES_HH
#define DUMUX_FV_GRID_VARIABLES_HH

#include <type_traits>
#include <memory>

// TODO: Remove once default template argument is omitted,
//       which is there solely to ensure backwards compatibility.
#include <dune/istl/bvector.hh>

#include <dumux/discretization/localview.hh>
#include "gridvariables.hh"

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Finite volume-specific local view on grid variables.
 * \tparam GV The grid variables class
 */
template<class GV>
class FVGridVariablesLocalView
{
    using GridGeometry = typename GV::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! export corresponding grid-wide class
    using GridVariables = GV;

    //! export underlying local views
    using ElementVolumeVariables = typename GV::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GV::GridFluxVariablesCache::LocalView;

    //! Constructor
    FVGridVariablesLocalView(const GridVariables& gridVariables)
    : gridVariables_(&gridVariables)
    , elemVolVars_(gridVariables.gridVolVars())
    , elemFluxVarsCache_(gridVariables.gridFluxVarsCache())
    {}

    /*!
     * \brief Bind this local view to a grid element.
     * \param element The grid element
     * \param fvGeometry Local view on the grid geometry
     */
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry)
    {
        const auto& x = gridVariables().dofs();
        elemVolVars_.bind(element, fvGeometry, x);
        elemFluxVarsCache_.bind(element, fvGeometry, elemVolVars_);
    }

    /*!
     * \brief Bind only the volume variables local view to a grid element.
     * \param element The grid element
     * \param fvGeometry Local view on the grid geometry
     */
    void bindElemVolVars(const Element& element,
                         const FVElementGeometry& fvGeometry)
    {
        elemVolVars_.bind(element, fvGeometry, gridVariables().dofs());

        // unbind flux variables cache
        elemFluxVarsCache_ = localView(gridVariables().gridFluxVarsCache());
    }

    //! return reference to the elem vol vars
    const ElementVolumeVariables& elemVolVars() const { return elemVolVars_; }
    ElementVolumeVariables& elemVolVars() { return elemVolVars_; }

    //! return reference to the flux variables cache
    const ElementFluxVariablesCache& elemFluxVarsCache() const { return elemFluxVarsCache_; }
    ElementFluxVariablesCache& elemFluxVarsCache() { return elemFluxVarsCache_; }

    //! Return reference to the grid variables
    const GridVariables& gridVariables() const
    { return *gridVariables_; }

private:
    const GridVariables* gridVariables_;
    ElementVolumeVariables elemVolVars_;
    ElementFluxVariablesCache elemFluxVarsCache_;
};

/*!
 * \ingroup Discretization
 * \brief The grid variable class for finite volume schemes, storing
 *        variables on scv and scvf (volume and flux variables).
 * \tparam GG the type of the grid geometry
 * \tparam GVV the type of the grid volume variables
 * \tparam GFVC the type of the grid flux variables cache
 * \tparam X the type used for solution vectors
 * \todo TODO: GG is an obsolote template argument?
 */
template<class GG, class GVV, class GFVC,
         class X = Dune::BlockVector<typename GVV::VolumeVariables::PrimaryVariables>>
class FVGridVariables
: public GridVariables<GG, X>
{
    using ParentType = GridVariables<GG, X>;
    using ThisType = FVGridVariables<GG, GVV, GFVC>;

public:
    using typename ParentType::SolutionVector;

    //! export type of the finite volume grid geometry
    using GridGeometry = GG;

    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    //! export primary variable type
    using PrimaryVariables = typename VolumeVariables::PrimaryVariables;

    //! export cache type for flux variables
    using GridFluxVariablesCache = GFVC;

    //! export the local view on this class
    using LocalView = FVGridVariablesLocalView<ThisType>;

    /*!
     * \brief Constructor
     * \param problem The problem to be solved
     * \param gridGeometry The geometry of the computational grid
     * \todo TODO: Here we could forward to the base class with
     *             [problem] (auto& x) { problem->applyInitialSolution(x); },
     *             but we cannot be sure that user problems implement
     *             the initial() or initialAtPos() functions?
     */
    template<class Problem>
    FVGridVariables(std::shared_ptr<Problem> problem,
                    std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , gridVolVars_(*problem)
    , prevGridVolVars_(*problem)
    , gridFluxVarsCache_(*problem)
    {}

    /*!
     * \brief Constructor with custom initialization of the solution.
     * \param problem The problem to be solved
     * \param gridGeometry The geometry of the computational grid
     * \param solOrInitializer This can be either a reference to a
     *                         solution vector, or an initializer
     *                         lambda. See Dumux::Variables.
     */
    template<class Problem, class SolOrInitializer>
    FVGridVariables(std::shared_ptr<Problem> problem,
                    std::shared_ptr<const GridGeometry> gridGeometry,
                    SolOrInitializer&& solOrInitializer)
    : ParentType(gridGeometry, std::forward<SolOrInitializer>(solOrInitializer))
    , gridVolVars_(*problem)
    , prevGridVolVars_(*problem)
    , gridFluxVarsCache_(*problem)
    {
        gridVolVars_.update(this->gridGeometry(), this->dofs());
        gridFluxVarsCache_.update(this->gridGeometry(), gridVolVars_, this->dofs(), true);
        prevGridVolVars_ = gridVolVars_;
    }

    //! initialize all variables (stationary case)
    [[deprecated("Initialize grid variables upon construction instead.")]]
    void init(const SolutionVector& curSol)
    {
        ParentType::update(curSol);

        // resize and update the volVars with the initial solution
        gridVolVars_.update(this->gridGeometry(), curSol);

        // update the flux variables caches (always force flux cache update on initialization)
        gridFluxVarsCache_.update(this->gridGeometry(), gridVolVars_, curSol, true);

        // set the volvars of the previous time step in case we have an instationary problem
        // note that this means some memory overhead in the case of enabled caching, however
        // this it outweighted by the advantage of having a single grid variables object for
        // stationary and instationary problems
        // TODO: Remove this in the context of new time integration schemes
        prevGridVolVars_ = gridVolVars_;
    }

    //! Deprecated update interface.
    [[deprecated("Use update(curSol) or forceUpdateAll(curSol) instead.")]]
    void update(const SolutionVector& curSol, bool forceFluxCacheUpdate)
    {
        if (forceFluxCacheUpdate)
            forceUpdateAll(curSol);
        else
            update(curSol);
    }

    //! update all variables after grid adaption
    [[deprecated("use forceUpdateAll() instead.")]]
    void updateAfterGridAdaption(const SolutionVector& curSol)
    { forceUpdateAll(curSol); }

    //! Update all variables that may be affected by a change in solution
    void update(const SolutionVector& curSol)
    {
        ParentType::update(curSol);

        // resize and update the volVars with the initial solution
        gridVolVars_.update(this->gridGeometry(), curSol);

        // update the flux variables caches
        gridFluxVarsCache_.update(this->gridGeometry(), gridVolVars_, curSol);
    }

    //! Force the update of all variables
    void forceUpdateAll(const SolutionVector& curSol)
    {
        ParentType::update(curSol);

        // resize and update the volVars with the initial solution
        gridVolVars_.update(this->gridGeometry(), curSol);

        // update the flux variables caches
        gridFluxVarsCache_.update(this->gridGeometry(), gridVolVars_, curSol, true);
    }

    /*!
     * \brief Sets the current state as the previous for next time step
     * \note this has to be called at the end of each time step
     */
    void advanceTimeStep()
    {
        prevGridVolVars_ = gridVolVars_;
    }

    //! resets state to the one before time integration
    void resetTimeStep(const SolutionVector& solution)
    {
        ParentType::update(solution);

        // set the new time step vol vars to old vol vars
        gridVolVars_ = prevGridVolVars_;

        // update the flux variables caches
        gridFluxVarsCache_.update(this->gridGeometry(), gridVolVars_, solution);
    }

    //! return the flux variables cache
    const GridFluxVariablesCache& gridFluxVarsCache() const
    { return gridFluxVarsCache_; }

    //! return the flux variables cache
    GridFluxVariablesCache& gridFluxVarsCache()
    { return gridFluxVarsCache_; }

    //! return the current volume variables
    const GridVolumeVariables& gridVolVars() const
    { return gridVolVars_; }

    //! return the current volume variables
    GridVolumeVariables& gridVolVars()
    { return gridVolVars_; }

    //! return the current volume variables
    [[deprecated("Use gridVolVars() instead")]]
    const GridVolumeVariables& curGridVolVars() const
    { return gridVolVars_; }

    //! return the current volume variables
    [[deprecated("Use gridVolVars() instead")]]
    GridVolumeVariables& curGridVolVars()
    { return gridVolVars_; }

    //! return the volume variables of the previous time step (for instationary problems)
    const GridVolumeVariables& prevGridVolVars() const
    { return prevGridVolVars_; }

    //! return the volume variables of the previous time step (for instationary problems)
    GridVolumeVariables& prevGridVolVars()
    { return prevGridVolVars_; }

private:
    GridVolumeVariables gridVolVars_;          //!< the current volume variables (primary and secondary variables)
    GridVolumeVariables prevGridVolVars_;      //!< the previous time step's volume variables (primary and secondary variables)
    GridFluxVariablesCache gridFluxVarsCache_; //!< the flux variables cache
};

} // end namespace Dumux

#endif
