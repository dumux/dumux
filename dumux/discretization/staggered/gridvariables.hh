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
#ifndef DUMUX_STAGGERED_GRID_VARIABLES_HH
#define DUMUX_STAGGERED_GRID_VARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/fvgridvariables.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Class storing scv and scvf variables
 */
template<class TypeTag>
class StaggeredGridVariables : public FVGridVariables<TypeTag>
{
    using ParentType = FVGridVariables<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using GridFaceVariables = typename GET_PROP_TYPE(TypeTag, GridFaceVariables);
    using GridFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

public:
    //! Constructor
    StaggeredGridVariables(std::shared_ptr<const Problem> problem,
                           std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(problem, fvGridGeometry)
    , fvGridGeometry_(fvGridGeometry)
    , curGridFaceVariables_(*problem)
    , prevGridFaceVariables_(*problem)
    {}

    //! update all variables
    void update(const SolutionVector& curSol)
    {
        ParentType::update(curSol);
        curGridFaceVariables_.update(*fvGridGeometry_, curSol);
    }

    //! initialize all variables (stationary case)
    void init(const SolutionVector& curSol)
    {
        ParentType::init(curSol);
        curGridFaceVariables_.update(*fvGridGeometry_, curSol);
    }

    //! initialize all variables (instationary case)
    void init(const SolutionVector& curSol, const SolutionVector& initSol)
    {
        ParentType::init(curSol, initSol);
        curGridFaceVariables_.update(*fvGridGeometry_, curSol);
        prevGridFaceVariables_.update(*fvGridGeometry_, initSol);
    }

    //! Sets the current state as the previous for next time step
    //! this has to be called at the end of each time step
    void advanceTimeStep()
    {
        ParentType::advanceTimeStep();
        prevGridFaceVariables_ = curGridFaceVariables_;
    }

    //! resets state to the one before time integration
    void resetTimeStep(const SolutionVector& solution)
    {
        ParentType::resetTimeStep(solution);
        curGridFaceVariables_ = prevGridFaceVariables_;
    }

    const GridFaceVariables& curGridFaceVars() const
    { return curGridFaceVariables_; }

    const GridFaceVariables& prevGridFaceVars() const
    { return prevGridFaceVariables_; }

    GridFaceVariables& curGridFaceVars()
    { return curGridFaceVariables_; }

    GridFaceVariables& prevGridFaceVars()
    { return prevGridFaceVariables_; }

private:

    std::shared_ptr<const FVGridGeometry> fvGridGeometry_;

    GridFaceVariables curGridFaceVariables_;
    GridFaceVariables prevGridFaceVariables_;
};

} // end namespace Dumux

#endif
