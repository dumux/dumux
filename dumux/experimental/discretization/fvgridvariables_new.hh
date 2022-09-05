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
#ifndef DUMUX_EXPERIMENTAL_DISCRETIZATION_FV_GRID_VARIABLES_HH
#define DUMUX_EXPERIMENTAL_DISCRETIZATION_FV_GRID_VARIABLES_HH

#include <utility>
#include <memory>
#include <type_traits>

#include <dumux/discretization/localview.hh>
#include <dumux/experimental/discretization/gridvariables.hh>

namespace Dumux::Experimental {

#ifndef DOXYGEN
namespace Detail {

template<typename P>
using ProblemGridGeometry = std::decay_t<decltype(std::declval<P>().gridGeometry())>;

} // namespace Detail
#endif // DOXYGEN

/*!
 * \ingroup Discretization
 * \brief Finite volume-specific local view on grid variables.
 * \tparam GV The grid variables class
 */
template<typename GV>
class FVGridVariablesLocalView
{
    using GridGeometry = typename GV::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using ElementVolumeVariables = typename GV::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GV::GridFluxVariablesCache::LocalView;

public:
    //! export corresponding grid-wide class
    using GridVariables = GV;

    //! Constructor
    FVGridVariablesLocalView(const GridVariables& gridVariables)
    : gridVariables_(&gridVariables)
    , elemVolVars_(gridVariables.gridVolVars())
    , elemFluxVarsCache_(gridVariables.gridFluxVarsCache())
    {}

    /*!
     * \brief Bind this local view to a grid element (r-value overload).
     * \param fvGeometry Local view on the grid geometry
     * \note Allows usage like `auto elemVars = localView(gridVars).bind(fvGeometry);`
     */
    FVGridVariablesLocalView bind(const FVElementGeometry& fvGeometry) &&
    {
        bind_(fvGeometry);
        return std::move(*this);
    }

    /*!
     * \brief Bind this local view to a grid element.
     * \param fvGeometry Local view on the grid geometry
     * \todo TODO: other bind policies!?
     */
    void bind(const FVElementGeometry& fvGeometry) &
    { bind_(fvGeometry); }

    const ElementVolumeVariables& elemVolVars() const { return elemVolVars_; }
    ElementVolumeVariables& elemVolVars() { return elemVolVars_; }

    const ElementFluxVariablesCache& elemFluxVarsCache() const { return elemFluxVarsCache_; }
    ElementFluxVariablesCache& elemFluxVarsCache() { return elemFluxVarsCache_; }

    const GridVariables& gridVariables() const
    { return *gridVariables_; }

private:
    void bind_(const FVElementGeometry& fvGeometry)
    {
        const auto& x = gridVariables().dofs();
        elemVolVars_.bind(fvGeometry.element(), fvGeometry, x);
        elemFluxVarsCache_.bind(fvGeometry.element(), fvGeometry, elemVolVars_);
    }

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
template<typename GVV,
         typename GFVC,
         typename X>
class FVGridVariables
: public GridVariables<Detail::ProblemGridGeometry<typename GVV::Problem>, X>
{
    using P = typename GVV::Problem;
    using GG = Detail::ProblemGridGeometry<typename GVV::Problem>;

    using ParentType = GridVariables<GG, X>;
    using ThisType = FVGridVariables<GVV, GFVC, X>;

public:
    using Problem = P;
    using GridGeometry = GG;
    using SolutionVector = X;
    using GridVolumeVariables = GVV;
    using GridFluxVariablesCache = GFVC;
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
    { updateCaches_(); }

    /*!
     * \brief Constructor with custom initialization of the solution.
     * \param problem The problem to be solved
     * \param gridGeometry The geometry of the computational grid
     * \param solOrInitializer This can be either a reference to a solution
     *                         vector, or an initializer lambda.
     *                         See Dumux::Experimental::Variables.
     */
    template<typename SolOrInitializer>
    FVGridVariables(std::shared_ptr<Problem> problem,
                    std::shared_ptr<const GridGeometry> gridGeometry,
                    SolOrInitializer&& solOrInitializer)
    : ParentType(gridGeometry, std::forward<SolOrInitializer>(solOrInitializer))
    , gridVolVars_(*problem)
    , gridFluxVarsCache_(*problem)
    { updateCaches_(); }

    //! Update the state to a new solution & time level
    template<typename Scalar>
    void update(const SolutionVector& x, const TimeLevel<Scalar>& t)
    {
        ParentType::updateTime(t);
        update(x);
    }

    //! Update all variables that may be affected by a change in solution
    void update(const SolutionVector& curSol)
    {
        ParentType::update(curSol);
        updateCaches_();
        // TODO: Invoke privarswitch
    }

    const GridFluxVariablesCache& gridFluxVarsCache() const
    { return gridFluxVarsCache_; }

    GridFluxVariablesCache& gridFluxVarsCache()
    { return gridFluxVarsCache_; }

    const GridVolumeVariables& gridVolVars() const
    { return gridVolVars_; }

    GridVolumeVariables& gridVolVars()
    { return gridVolVars_; }

private:
    void updateCaches_()
    {
        gridVolVars_.update(this->gridGeometry(), this->dofs());
        gridFluxVarsCache_.update(this->gridGeometry(), gridVolVars_, this->dofs());
    }

    GridVolumeVariables gridVolVars_;          //!< the current volume variables (primary and secondary variables)
    GridFluxVariablesCache gridFluxVarsCache_; //!< the flux variables cache
};

} // end namespace Dumux::Experimental

#endif
