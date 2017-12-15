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
 * \brief Class storing scv and scvf variables
 */
#ifndef DUMUX_FV_GRID_VARIABLES_HH
#define DUMUX_FV_GRID_VARIABLES_HH

#include <memory>
#include <dumux/common/properties.hh>

namespace Dumux
{

/*!
 * \brief Class storing scv and scvf variables
 */
template<class TypeTag>
class FVGridVariables
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using GridFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GlobalFluxVariablesCache);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

public:
    //! Constructor
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

    //! Sets the current state as the previous for next time step
    //! this has to be called at the end of each time step
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

    const GridFluxVariablesCache& gridFluxVarsCache() const
    { return gridFluxVarsCache_; }

    const GridVolumeVariables& curGridVolVars() const
    { return curGridVolVars_; }

    const GridVolumeVariables& prevGridVolVars() const
    { return prevGridVolVars_; }

    GridFluxVariablesCache& gridFluxVarsCache()
    { return gridFluxVarsCache_; }

    GridVolumeVariables& curGridVolVars()
    { return curGridVolVars_; }

    GridVolumeVariables& prevGridVolVars()
    { return prevGridVolVars_; }

private:
    std::shared_ptr<const Problem> problem_;
    std::shared_ptr<const FVGridGeometry> fvGridGeometry_;

    // the current and previous variables (primary and secondary variables)
    GridVolumeVariables curGridVolVars_;
    GridVolumeVariables prevGridVolVars_;

    // the flux variables cache vector vector
    GridFluxVariablesCache gridFluxVarsCache_;
};

} // end namespace Dumux

#endif
