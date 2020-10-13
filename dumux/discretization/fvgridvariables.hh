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
 * \brief The grid variable class for finite volume schemes storing variables on scv and scvf (volume and flux variables)
 */
#ifndef DUMUX_FV_GRID_VARIABLES_HH
#define DUMUX_FV_GRID_VARIABLES_HH

#include <type_traits>
#include <memory>
#include <cassert>

#include <dumux/timestepping/timelevel.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Local view on the grid variables for finite volume schemes.
 * \note This default implementation only stores a pointer on the grid variables.
 * \tparam GV The grid variables class
 */
template<class GV>
class FVGridVariablesLocalView
{
    using GridGeometry = typename GV::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using ElementVolumeVariables = typename GV::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GV::GridFluxVariablesCache::LocalView;

public:
    using GridVariables = GV;

    //! The constructor
    FVGridVariablesLocalView(const GridVariables& gridVariables)
    : gridVariables_(&gridVariables)
    , elemVolVars_(gridVariables.gridVolVars())
    , elemFluxVarsCache_(gridVariables.gridFluxVarsCache())
    {}

    /*!
     * \brief Bind this local view to a grid element.
     * \param element The grid element
     * \param fvGeometry Local view on the grid geometry
     * \todo TODO: We probably need selective binds also, as sometimes
     *             one does not want to bind the flux cache, but only the
     *             vol vars? For instance, when computing masses or storage
     *             after a time step...
     */
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry)
    {
        const auto& x = gridVariables().dofs();

        elemVolVars_.bind(element, fvGeometry, x);
        elemFluxVarsCache_.bind(element, fvGeometry, elemVolVars_);
    }

    // TODO: Doc me
    const ElementVolumeVariables& elemVolVars() const { return elemVolVars_; }
    const ElementFluxVariablesCache& elemFluxVarsCache() const { return elemFluxVarsCache_; }

    // TODO: Doc me
    // We need those non-const variants to be able to deflect things in numeric assembly
    ElementVolumeVariables& elemVolVars() { return elemVolVars_; }
    ElementFluxVariablesCache& elemFluxVarsCache() { return elemFluxVarsCache_; }


    /*!
     * \brief Return a reference to the grid variables.
     */
    const GridVariables& gridVariables() const
    { return *gridVariables_; }

private:
    const GridVariables* gridVariables_;
    ElementVolumeVariables elemVolVars_;
    ElementFluxVariablesCache elemFluxVarsCache_;
};

/*!
 * \ingroup Discretization
 * \brief The grid variable class for finite volume schemes storing variables on scv and scvf (volume and flux variables)
 * \tparam the type of the grid geometry
 * \tparam the type of the grid volume variables
 * \tparam the type of the grid flux variables cache
 * \todo TODO: Should Problem be template arg instead of GG? We expect a problem
 *             to be defined on a grid geometry anyway, so they always carry the
 *             grid geometry with them...
 * \todo TODO: We also need to know the SolutionVector type...
 */
template<class GG, class GVV, class GFVC>
class FVGridVariables
{
    using ThisType = FVGridVariables<GG, GVV, GFVC>;

public:
    //! export type of the finite volume grid geometry
    using GridGeometry = GG;

    //! export type of the finite volume grid geometry
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    //! export primary variable type
    //! TODO: Must this always match the value_type of solution vector?
    //!       In that case we should check maybe...
    using PrimaryVariables = typename VolumeVariables::PrimaryVariables;

    //! export scalar type (TODO get it directly from the volvars)
    using Scalar = std::decay_t<decltype(std::declval<PrimaryVariables>()[0])>;

    //! export type of the finite volume grid geometry
    using GridFluxVariablesCache = GFVC;

    //! The local view on these grid variables
    using LocalView = FVGridVariablesLocalView<ThisType>;

    //! export the underlying problem to be solved
    using Problem = typename GVV::Problem;

    //! TODO: SolutionVector as template arg?
    using SolutionVector = Dune::BlockVector<PrimaryVariables>;

    //! Type representing a time level
    using TimeLevel = Dumux::TimeLevel<Scalar>;

    //! constructor initializing the solution from the problem
    // TODO: Constructor with lambda (see FEGridVariables)
    template<class Problem>
    FVGridVariables(std::shared_ptr<Problem> problem,
                    const TimeLevel& timeLevel = TimeLevel{0.0})
    : problem_(problem)
    , gridVolVars_(*problem)
    , gridFluxVarsCache_(*problem)
    , timeLevel_(timeLevel)
    {
        problem->applyInitialSolution(x_);
        gridVolVars_.update(problem->gridGeometry(), x_);
        gridFluxVarsCache_.update(problem->gridGeometry(), gridVolVars_, x_, /*forceUpdate???*/true);
    }

    //! return the underlying problem
    const Problem& problem() const
    { return *problem_; }

    //! return the finite volume grid geometry
    const GridGeometry& gridGeometry() const
    { return problem_->gridGeometry(); }

    //! return the solution for which the grid variables were updated
    const SolutionVector& dofs() const
    { return x_; }

    //! return the time level the grid variables represent
    const TimeLevel& timeLevel() const
    { return timeLevel_; }

    //! return the flux variables cache
    const GridFluxVariablesCache& gridFluxVarsCache() const
    { return gridFluxVarsCache_; }

    //! return the flux variables cache
    //! TODO: This is needed for the case of caching! Can we circumvent this???
    GridFluxVariablesCache& gridFluxVarsCache()
    { return gridFluxVarsCache_; }

    //! return the current volume variables
    const GridVolumeVariables& gridVolVars() const
    { return gridVolVars_; }

    //! return the current volume variables
    //! TODO: This is needed for the case of caching! Can we circumvent this???
    GridVolumeVariables& gridVolVars()
    { return gridVolVars_; }

    //! update all variables subject to a new solution
    void updateDofs(const SolutionVector& x)
    {
        x_ = x;
        gridVolVars_.update(problem().gridGeometry(), x_);
        gridFluxVarsCache_.update(problem().gridGeometry(), gridVolVars_, x_);
    }

    //! TODO Doc me
    void updateTime(const TimeLevel& timeLevel)
    {
        timeLevel_ = timeLevel;
    }

    //! update all variables subject to a new solution
    void update(const SolutionVector& x, const TimeLevel& timeLevel)
    {
        updateDofs(x);
        updateTime(timeLevel);
    }

private:
    std::shared_ptr<const Problem> problem_; //!< pointer to the problem to be solved

    GridVolumeVariables gridVolVars_;          //!< the volume variables (primary and secondary variables)
    GridFluxVariablesCache gridFluxVarsCache_; //!< the flux variables cache

    SolutionVector x_;    //!< the solution corresponding to this class' current state
    TimeLevel timeLevel_; //!< info during time integration for instationary problems
};

} // end namespace Dumux

#endif
