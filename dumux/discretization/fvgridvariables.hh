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
 * \ingroup Discretization
 * \brief The grid variable class for finite volume schemes storing variables on scv and scvf (volume and flux variables)
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
                    std::shared_ptr<GridGeometry> fvGridGeometry)
    : fvGridGeometry_(fvGridGeometry)
    , curGridVolVars_(*problem)
    , prevGridVolVars_(*problem)
    , gridFluxVarsCache_(*problem)
    {}

    //! update all variables
    template<class SolutionVector>
    void update(const SolutionVector& curSol, bool forceFluxCacheUpdate = false)
    {
        // resize and update the volVars with the initial solution
        curGridVolVars_.update(*fvGridGeometry_, curSol);

        // update the flux variables caches
        gridFluxVarsCache_.update(*fvGridGeometry_, curGridVolVars_, curSol, forceFluxCacheUpdate);
    }

    //! update all variables after grid adaption
    template<class SolutionVector>
    void updateAfterGridAdaption(const SolutionVector& curSol)
    {
        // update (always force flux cache update as the grid changed)
        update(curSol, true);

        // for instationary problems also update the variables
        // for the previous time step to the new grid
        if (!problemIsStationary_)
            prevGridVolVars_ = curGridVolVars_;
    }

    //! initialize all variables (stationary case)
    template<class SolutionVector>
    void init(const SolutionVector& curSol)
    {
        // resize and update the volVars with the initial solution
        curGridVolVars_.update(*fvGridGeometry_, curSol);

        // update the flux variables caches (always force flux cache update on initialization)
        gridFluxVarsCache_.update(*fvGridGeometry_, curGridVolVars_, curSol, true);
    }

    //! initialize all variables (instationary case)
    template<class SolutionVector>
    void init(const SolutionVector& curSol, const SolutionVector& initSol)
    {
        // remember that we have a stationary problem
        problemIsStationary_ = false;

        // initialize current volvars and the flux var cache
        init(curSol);

        // update the old time step vol vars with the initial solution
        prevGridVolVars_.update(*fvGridGeometry_, initSol);
    }

    /*!
     * \brief Sets the current state as the previous for next time step
     * \note this has to be called at the end of each time step
     */
    void advanceTimeStep()
    {
        assert(!problemIsStationary_);
        prevGridVolVars_ = curGridVolVars_;
    }

    //! resets state to the one before time integration
    template<class SolutionVector>
    void resetTimeStep(const SolutionVector& solution)
    {
        assert(!problemIsStationary_);

        // set the new time step vol vars to old vol vars
        curGridVolVars_ = prevGridVolVars_;

        // update the flux variables caches
        gridFluxVarsCache_.update(*fvGridGeometry_, curGridVolVars_, solution);
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

protected:

    std::shared_ptr<const GridGeometry> fvGridGeometry_; //!< pointer to the constant grid geometry

private:
    GridVolumeVariables curGridVolVars_; //!< the current volume variables (primary and secondary variables)
    GridVolumeVariables prevGridVolVars_; //!< the previous time step's volume variables (primary and secondary variables)

    GridFluxVariablesCache gridFluxVarsCache_; //!< the flux variables cache

    bool problemIsStationary_ = true;
};

} // end namespace Dumux

#endif
