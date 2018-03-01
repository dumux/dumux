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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredGridVariables
 */
#ifndef DUMUX_STAGGERED_GRID_VARIABLES_HH
#define DUMUX_STAGGERED_GRID_VARIABLES_HH

#include <dumux/discretization/fvgridvariables.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class storing data associated to scvs and scvfs
 * \tparam GG the type of the grid geometry
 * \tparam GVV the type of the grid volume variables
 * \tparam GFVC the type of the grid flux variables cache
 * \tparam GFV the type of the grid face variables
 */
template<class GG, class GVV, class GFVC, class GFV>
class StaggeredGridVariables : public FVGridVariables<GG, GVV, GFVC>
{
    using ParentType = FVGridVariables<GG, GVV, GFVC>;
    using FVGridGeometry = GG;
    using GridVolumeVariables = GVV;
    using GridFluxVariablesCache = GFVC;
    using GridFaceVariables = GFV;

public:
    //! Constructor
    template<class Problem>
    StaggeredGridVariables(std::shared_ptr<Problem> problem,
                           std::shared_ptr<FVGridGeometry> fvGridGeometry)
    : ParentType(problem, fvGridGeometry)
    , curGridFaceVariables_(*problem)
    , prevGridFaceVariables_(*problem)
    {}

    //! update all variables
    template<class SolutionVector>
    void update(const SolutionVector& curSol)
    {
        ParentType::update(curSol);
        curGridFaceVariables_.update(*this->fvGridGeometry_, curSol);
    }

    //! initialize all variables (stationary case)
    template<class SolutionVector>
    void init(const SolutionVector& curSol)
    {
        ParentType::init(curSol);
        curGridFaceVariables_.update(*this->fvGridGeometry_, curSol);
    }

    //! initialize all variables (instationary case)
    template<class SolutionVector>
    void init(const SolutionVector& curSol, const SolutionVector& initSol)
    {
        ParentType::init(curSol, initSol);
        curGridFaceVariables_.update(*this->fvGridGeometry_, curSol);
        prevGridFaceVariables_.update(*this->fvGridGeometry_, initSol);
    }

    //! Sets the current state as the previous for next time step
    //! this has to be called at the end of each time step
    void advanceTimeStep()
    {
        ParentType::advanceTimeStep();
        prevGridFaceVariables_ = curGridFaceVariables_;
    }

    //! resets state to the one before time integration
    template<class SolutionVector>
    void resetTimeStep(const SolutionVector& solution)
    {
        ParentType::resetTimeStep(solution);
        curGridFaceVariables_ = prevGridFaceVariables_;
    }

    //! return the current face variables
    const GridFaceVariables& curGridFaceVars() const
    { return curGridFaceVariables_; }

    //! return the previous face variables
    const GridFaceVariables& prevGridFaceVars() const
    { return prevGridFaceVariables_; }

    //! return the current face variables
    GridFaceVariables& curGridFaceVars()
    { return curGridFaceVariables_; }

    //! return the previous face variables
    GridFaceVariables& prevGridFaceVars()
    { return prevGridFaceVariables_; }

private:
    GridFaceVariables curGridFaceVariables_;
    GridFaceVariables prevGridFaceVariables_;
};

} // end namespace Dumux

#endif
