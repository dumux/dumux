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

#include <memory>
#include <dumux/common/properties.hh>

namespace Dumux
{

/*!
 * \ingroup Discretization
 * \brief The grid variable class for finite volume schemes storing variables on scv and scvf (volume and flux variables)
 */
template<class TypeTag>
class FVGridVariables
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using GridFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

public:
    FVGridVariables(std::shared_ptr<const Problem> problem,
                    std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : problem_(problem)
    , fvGridGeometry_(fvGridGeometry)
    , curGridVolVars_(*problem)
    , prevGridVolVars_(*problem)
    , gridFluxVarsCache_(*problem)
    {}

    //! update all variables
    void update(const SolutionVector& curSol)
    {
        // resize and update the volVars with the initial solution
        curGridVolVars_.update(*fvGridGeometry_, curSol);

        // update the flux variables caches
        gridFluxVarsCache_.update(*fvGridGeometry_, curGridVolVars_, curSol);
    }

    //! initialize all variables (stationary case)
    void init(const SolutionVector& curSol)
    {
        // resize and update the volVars with the initial solution
        curGridVolVars_.update(*fvGridGeometry_, curSol);

        // update the flux variables caches
        gridFluxVarsCache_.update(*fvGridGeometry_, curGridVolVars_, curSol, true);
    }

    //! initialize all variables (instationary case)
    void init(const SolutionVector& curSol, const SolutionVector& initSol)
    {
        // resize and update the volVars with the initial solution
        curGridVolVars_.update(*fvGridGeometry_, curSol);

        // update the flux variables caches
        gridFluxVarsCache_.update(*fvGridGeometry_, curGridVolVars_, curSol, true);

        // update the old time step vol vars with the initial solution
        prevGridVolVars_.update(*fvGridGeometry_, initSol);
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
    void resetTimeStep(const SolutionVector& solution)
    {
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

private:
    std::shared_ptr<const Problem> problem_; //!< pointer to the constant problem definition
    std::shared_ptr<const FVGridGeometry> fvGridGeometry_; //!< pointer to the constant grid geometry

    GridVolumeVariables curGridVolVars_; //!< the current volume variables (primary and secondary variables)
    GridVolumeVariables prevGridVolVars_; //!< the previous time step's volume variables (primary and secondary variables)

    GridFluxVariablesCache gridFluxVarsCache_; //!< the flux variables cache
};

} // end namespace Dumux

#endif
