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

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief The grid variable class for finite volume schemes storing variables on scv and scvf (volume and flux variables)
 * \tparam the type of the grid geometry
 * \tparam the type of the grid volume variables
 * \tparam the type of the grid flux variables cache
 */
template<class GG, class GVV, class GFVC>
class FVGridVariables
{
public:
    //! export type of the finite volume grid geometry
    using GridGeometry = GG;

    //! export type of the finite volume grid geometry
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    //! export primary variable type
    using PrimaryVariables = typename VolumeVariables::PrimaryVariables;

    //! export scalar type (TODO get it directly from the volvars)
    using Scalar = std::decay_t<decltype(std::declval<PrimaryVariables>()[0])>;

    //! export type of the finite volume grid geometry
    using GridFluxVariablesCache = GFVC;

    template<class Problem>
    FVGridVariables(std::shared_ptr<Problem> problem,
                    std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_(gridGeometry)
    , curGridVolVars_(*problem)
    , prevGridVolVars_(*problem)
    , gridFluxVarsCache_(*problem)
    {}

    //! initialize all variables (stationary case)
    template<class SolutionVector>
    void init(const SolutionVector& curSol)
    {
        // resize and update the volVars with the initial solution
        curGridVolVars_.update(*gridGeometry_, curSol);

        // update the flux variables caches (always force flux cache update on initialization)
        gridFluxVarsCache_.update(*gridGeometry_, curGridVolVars_, curSol, true);

        // set the volvars of the previous time step in case we have an instationary problem
        // note that this means some memory overhead in the case of enabled caching, however
        // this it outweighted by the advantage of having a single grid variables object for
        // stationary and instationary problems
        prevGridVolVars_ = curGridVolVars_;
    }

    //! update all variables
    template<class SolutionVector>
    void update(const SolutionVector& curSol, bool forceFluxCacheUpdate = false)
    {
        // resize and update the volVars with the initial solution
        curGridVolVars_.update(*gridGeometry_, curSol);

        // update the flux variables caches
        gridFluxVarsCache_.update(*gridGeometry_, curGridVolVars_, curSol, forceFluxCacheUpdate);
    }

    //! update all variables after grid adaption
    template<class SolutionVector>
    void updateAfterGridAdaption(const SolutionVector& curSol)
    {
        // update (always force flux cache update as the grid changed)
        update(curSol, true);

        // for instationary problems also update the variables
        // for the previous time step to the new grid
        prevGridVolVars_ = curGridVolVars_;
    }

    /*!
     * \brief Sets the current state as the previous for next time step
     * \note this has to be called at the end of each time step
     */
    void advanceTimeStep()
    {
        prevGridVolVars_ = curGridVolVars_;
    }

    //! resets state to the one before time integration
    template<class SolutionVector>
    void resetTimeStep(const SolutionVector& solution)
    {
        // set the new time step vol vars to old vol vars
        curGridVolVars_ = prevGridVolVars_;

        // update the flux variables caches
        gridFluxVarsCache_.update(*gridGeometry_, curGridVolVars_, solution);
    }

    //! return the flux variables cache
    const GridFluxVariablesCache& gridFluxVarsCache() const
    { return gridFluxVarsCache_; }

    //! return the flux variables cache
    GridFluxVariablesCache& gridFluxVarsCache()
    { return gridFluxVarsCache_; }

    //! return the current volume variables
    const GridVolumeVariables& curGridVolVars() const
    { return curGridVolVars_; }

    //! return the current volume variables
    GridVolumeVariables& curGridVolVars()
    { return curGridVolVars_; }

    //! return the volume variables of the previous time step (for instationary problems)
    const GridVolumeVariables& prevGridVolVars() const
    { return prevGridVolVars_; }

    //! return the volume variables of the previous time step (for instationary problems)
    GridVolumeVariables& prevGridVolVars()
    { return prevGridVolVars_; }

    //! return the finite volume grid geometry
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

protected:

    std::shared_ptr<const GridGeometry> gridGeometry_; //!< pointer to the constant grid geometry

private:
    GridVolumeVariables curGridVolVars_; //!< the current volume variables (primary and secondary variables)
    GridVolumeVariables prevGridVolVars_; //!< the previous time step's volume variables (primary and secondary variables)

    GridFluxVariablesCache gridFluxVarsCache_; //!< the flux variables cache
};

} // end namespace Dumux

//////////////////////////////////////////////////////////////
// Experimental implementation of new grid variables layout //
//////////////////////////////////////////////////////////////

#include <utility>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/discretization/localview.hh>
#include <dumux/discretization/gridvariables.hh>
#include <dumux/discretization/solutionstate.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Discretization
 * \brief Policy for binding local views of finite volume grid variables.
 */
enum class FVGridVariablesBindPolicy
{
    all,
    volVarsOnly,
    volVarsOfElementOnly
};

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

    //! export underlying volume variables cache
    using ElementVolumeVariables = typename GV::GridVolumeVariables::LocalView;

    //! export underlying flux variables cache
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
              const FVElementGeometry& fvGeometry,
              FVGridVariablesBindPolicy policy = FVGridVariablesBindPolicy::all)
    {
        const auto solState = solutionStateView(gridVariables());

        if (policy == FVGridVariablesBindPolicy::all)
        {
            elemVolVars_.bind(element, fvGeometry, solState);
            elemFluxVarsCache_.bind(element, fvGeometry, elemVolVars_);
        }
        else
        {
            // unbind flux variables cache
            elemFluxVarsCache_ = localView(gridVariables().gridFluxVarsCache());
            if (policy == FVGridVariablesBindPolicy::volVarsOnly)
                elemVolVars_.bind(element, fvGeometry, solState);
            else if (policy == FVGridVariablesBindPolicy::volVarsOfElementOnly)
                elemVolVars_.bindElement(element, fvGeometry, solState);
            else
                DUNE_THROW(Dune::NotImplemented, "Bind for given policy");
        }
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
 * \tparam GVV the type of the grid volume variables
 * \tparam GFVC the type of the grid flux variables cache
 * \tparam X the type used for solution vectors
 */
template<class GVV, class GFVC, class X>
class FVGridVariables
: public GridVariables<typename ProblemTraits<typename GVV::Problem>::GridGeometry, X>
{
    using Problem = typename GVV::Problem;
    using GG = typename ProblemTraits<Problem>::GridGeometry;

    using ParentType = GridVariables<GG, X>;
    using ThisType = FVGridVariables<GVV, GFVC, X>;

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
     * \note This constructor initializes the solution using the
     *       initializer function in the given problem, and thus,
     *       this only compiles if the problem implements it.
     */
    FVGridVariables(std::shared_ptr<Problem> problem,
                    std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry, [problem] (auto& x) { problem->applyInitialSolution(x); })
    , gridVolVars_(*problem)
    , gridFluxVarsCache_(*problem)
    {}

    /*!
     * \brief Constructor with custom initialization of the solution.
     * \param problem The problem to be solved
     * \param gridGeometry The geometry of the computational grid
     * \param solOrInitializer This can be either a reference to a solution
     *                         vector, or an initializer lambda.
     *                         See Dumux::Experimental::Variables.
     */
    template<class SolOrInitializer>
    FVGridVariables(std::shared_ptr<Problem> problem,
                    std::shared_ptr<const GridGeometry> gridGeometry,
                    SolOrInitializer&& solOrInitializer)
    : ParentType(gridGeometry, std::forward<SolOrInitializer>(solOrInitializer))
    , gridVolVars_(*problem)
    , gridFluxVarsCache_(*problem)
    {
        gridVolVars_.update(this->gridGeometry(), this->dofs());
        gridFluxVarsCache_.update(this->gridGeometry(), gridVolVars_, this->dofs(), true);
    }

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

private:
    GridVolumeVariables gridVolVars_;          //!< the current volume variables (primary and secondary variables)
    GridFluxVariablesCache gridFluxVarsCache_; //!< the flux variables cache
};

    // implementation details
    namespace Detail {
        struct hasGridVolVars
        {
            template<class GV>
            auto operator()(const GV& gv) -> decltype(gv.gridVolVars()) {}
        };
    } // end namespace Detail

    // helper bool to check if a grid variables type fulfills the new experimental interface
    // can be used in the transition period in places where both interfaces should be supported
    // Note that this doesn't check if the provided type actually models grid variables. Instead,
    // it is assumed that that is known, and it is only checked if the new or old behaviour can be
    // expected from it.
    template<class GridVariables>
    inline constexpr bool areExperimentalGridVars =
        decltype(isValid(Detail::hasGridVolVars())(std::declval<GridVariables>()))::value;

} // end namespace Dumux::Experimental

#endif
