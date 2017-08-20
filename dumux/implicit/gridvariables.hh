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
#ifndef DUMUX_GRID_VARIABLES_HH
#define DUMUX_GRID_VARIABLES_HH

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Class storing scv and scvf variables
 */
template<class TypeTag>
class GridVariables
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GlobalVolumeVariables);
    using GridFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GlobalFluxVariablesCache);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

public:
    //! Constructor
    GridVariables() {}

    //! update all fvElementGeometries (do this again after grid adaption)
    void update(const Problem& problem,
                const FVGridGeometry& fvGridGeometry,
                const SolutionVector& solution)
    {
        // resize and update the volVars with the initial solution
        curGlobalVolVars_.update(problem, fvGridGeometry, solution);

        // update the flux variables caches
        globalfluxVarsCache_.update(problem, fvGridGeometry, curGlobalVolVars_, solution);

        // set the old time step vol vars to new current vol vars
        prevGlobalVolVars_ = curGlobalVolVars_;
    }

    const GridFluxVariablesCache& gridFluxVarsCache() const
    { return globalfluxVarsCache_; }

    const GridVolumeVariables& curGridVolVars() const
    { return curGlobalVolVars_; }

    const GridVolumeVariables& prevGridVolVars() const
    { return prevGlobalVolVars_; }

    GridFluxVariablesCache& gridFluxVarsCache()
    { return globalfluxVarsCache_; }

    GridVolumeVariables& curGridVolVars()
    { return curGlobalVolVars_; }

    GridVolumeVariables& prevGridVolVars()
    { return prevGlobalVolVars_; }

private:

    // the current and previous variables (primary and secondary variables)
    GridVolumeVariables curGlobalVolVars_;
    GridVolumeVariables prevGlobalVolVars_;

    // the flux variables cache vector vector
    GridFluxVariablesCache globalfluxVarsCache_;
};

} // end namespace

#endif
