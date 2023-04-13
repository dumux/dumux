// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \ingroup Discretization
 * \brief The grid variable class for finite volume schemes,
 *        storing variables on scv and scvf (volume and flux variables)
 */
#ifndef DUMUX_EXPERIMENTAL_FV_GRID_VARIABLES_HH
#define DUMUX_EXPERIMENTAL_FV_GRID_VARIABLES_HH

#include <utility>
#include <memory>

#include <dumux/common/typetraits/problem.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/experimental/discretization/gridvariables.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
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
 * \ingroup Experimental
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

} // end namespace Dumux::Experimental

#endif
