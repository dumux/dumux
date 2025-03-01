// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief The grid variable class for cvfe schemes,
 *        storing variables on local dofs and scvf
 */
#ifndef DUMUX_CVFE_GRID_VARIABLES_HH
#define DUMUX_CVFE_GRID_VARIABLES_HH

#include <type_traits>
#include <memory>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief The grid variable class for cvfe schemes storing variables on local dofs and scvf
 * \tparam the type of the grid geometry
 * \tparam the type of the grid local variables
 * \tparam the type of the grid flux variables cache
 */
template<class GG, class GLV, class GFVC>
class CVFEGridVariables
{
public:
    //! export type of the grid geometry
    using GridGeometry = GG;

    //! export type of the grid local variables
    using GridLocalVariables = GLV;

    //! export the grid volume variables type
    using GridVolumeVariables [[deprecated("Use GridLocalVariables instead. Will be removed after release 3.10.")]]= GridLocalVariables;

    //! export type of the local variables
    using LocalVariables = typename GridLocalVariables::LocalVariables;

    //! export the volume variables type
    using VolumeVariables [[deprecated("Use LocalVariables instead. Will be removed after release 3.10.")]] = LocalVariables;

    //! export primary variable type
    using PrimaryVariables = typename LocalVariables::PrimaryVariables;

    //! export scalar type (TODO get it directly from the volvars)
    using Scalar = std::decay_t<decltype(std::declval<PrimaryVariables>()[0])>;

    //! export type of the cache for grid flux variables
    using GridFluxVariablesCache = GFVC;

    template<class Problem>
    CVFEGridVariables(std::shared_ptr<Problem> problem,
                    std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_(gridGeometry)
    , curGridLocalVars_(*problem)
    , prevGridLocalVars_(*problem)
    , gridFluxVarsCache_(*problem)
    {}

    //! initialize all variables (stationary case)
    template<class SolutionVector>
    void init(const SolutionVector& curSol)
    {
        // resize and update the volVars with the initial solution
        curGridLocalVars_.update(*gridGeometry_, curSol);

        // update the flux variables caches (always force flux cache update on initialization)
        gridFluxVarsCache_.update(*gridGeometry_, curGridLocalVars_, curSol, true);

        // set the volvars of the previous time step in case we have an instationary problem
        // note that this means some memory overhead in the case of enabled caching, however
        // this it outweighted by the advantage of having a single grid variables object for
        // stationary and instationary problems
        prevGridLocalVars_ = curGridLocalVars_;
    }

    //! update all variables
    template<class SolutionVector>
    void update(const SolutionVector& curSol, bool forceFluxCacheUpdate = false)
    {
        // resize and update the volVars with the initial solution
        curGridLocalVars_.update(*gridGeometry_, curSol);

        // update the flux variables caches
        gridFluxVarsCache_.update(*gridGeometry_, curGridLocalVars_, curSol, forceFluxCacheUpdate);
    }

    //! update all variables after grid adaption
    template<class SolutionVector>
    void updateAfterGridAdaption(const SolutionVector& curSol)
    {
        // update (always force flux cache update as the grid changed)
        update(curSol, true);

        // for instationary problems also update the variables
        // for the previous time step to the new grid
        prevGridLocalVars_ = curGridLocalVars_;
    }

    /*!
     * \brief Sets the current state as the previous for next time step
     * \note this has to be called at the end of each time step
     */
    void advanceTimeStep()
    {
        prevGridLocalVars_ = curGridLocalVars_;
    }

    //! resets state to the one before time integration
    template<class SolutionVector>
    void resetTimeStep(const SolutionVector& solution)
    {
        // set the new time step vol vars to old vol vars
        curGridLocalVars_ = prevGridLocalVars_;

        // update the flux variables caches
        gridFluxVarsCache_.update(*gridGeometry_, curGridLocalVars_, solution);
    }

    //! return the flux variables cache
    const GridFluxVariablesCache& gridFluxVarsCache() const
    { return gridFluxVarsCache_; }

    //! return the flux variables cache
    GridFluxVariablesCache& gridFluxVarsCache()
    { return gridFluxVarsCache_; }

    //! return the current local variables
    const GridLocalVariables& curGridLocalVars() const
    { return curGridLocalVars_; }

    //! return the current local variables
    GridLocalVariables& curGridLocalVars()
    { return curGridLocalVars_; }

    //! return the local variables of the previous time step (for instationary problems)
    const GridLocalVariables& prevGridLocalVars() const
    { return prevGridLocalVars_; }

    //! return the local variables of the previous time step (for instationary problems)
    GridLocalVariables& prevGridLocalVars()
    { return prevGridLocalVars_; }

    //! return the current local variables
    [[deprecated("Use curGridLocalVars instead. Will be removed after release 3.10.")]]
    const GridLocalVariables& curGridVolVars() const
    { return curGridLocalVars_; }

    //! return the current local variables
    [[deprecated("Use curGridLocalVars instead. Will be removed after release 3.10.")]]
    GridLocalVariables& curGridVolVars()
    { return curGridLocalVars_; }

    //! return the local variables of the previous time step (for instationary problems)
    [[deprecated("Use prevGridLocalVars instead. Will be removed after release 3.10.")]]
    const GridLocalVariables& prevGridVolVars() const
    { return prevGridLocalVars_; }

    //! return the local variables of the previous time step (for instationary problems)
    [[deprecated("Use prevGridLocalVars instead. Will be removed after release 3.10.")]]
    GridLocalVariables& prevGridVolVars()
    { return prevGridLocalVars_; }

    //! return the finite local grid geometry
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

protected:

    std::shared_ptr<const GridGeometry> gridGeometry_; //!< pointer to the constant grid geometry

private:
    GridLocalVariables curGridLocalVars_; //!< the current local variables (primary and secondary variables)
    GridLocalVariables prevGridLocalVars_; //!< the previous time step's local variables (primary and secondary variables)

    GridFluxVariablesCache gridFluxVarsCache_; //!< the flux variables cache
};

} // end namespace Dumux

#endif
