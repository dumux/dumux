// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief The grid variables class for for general schemes,
 *        storing variables and data on the grid
 */
#ifndef DUMUX_GRID_VARIABLES_HH
#define DUMUX_GRID_VARIABLES_HH

#include <type_traits>
#include <memory>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief The grid variables class for general schemes, storing variables and data on the grid
 * \tparam the type of the grid geometry
 * \tparam the type of the grid variables cache
 */
template<class GG, class GVC>
class GridVariables
{
public:
    //! export type of the finite volume grid geometry
    using GridGeometry = GG;

    //! export type of the grid variables cache
    using GridVariablesCache = GVC;

    //! export type of the volume variables
    using Variables = typename GridVariablesCache::Variables;

    //! export primary variable type
    using PrimaryVariables = typename Variables::PrimaryVariables;

    //! export scalar
    using Scalar = typename PrimaryVariables::value_type;

    template<class Problem>
    GridVariables(std::shared_ptr<Problem> problem,
                  std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_(gridGeometry)
    , curGridVars_(*problem)
    , prevGridVars_(*problem)
    {}

    //! initialize all variables (stationary case)
    template<class SolutionVector>
    void init(const SolutionVector& curSol)
    {
        // resize and update the volVars with the initial solution
        curGridVars_.init(*gridGeometry_, curSol);

        // set the variables of the previous time step in case we have an instationary problem
        // note that this means some memory overhead in the case of enabled caching, however
        // this it outweighted by the advantage of having a single grid variables object for
        // stationary and instationary problems
        prevGridVars_ = curGridVars_;
    }

    //! update all variables
    template<class SolutionVector>
    void update(const SolutionVector& curSol)
    {
        // resize and update the volVars with the initial solution
        curGridVars_.update(*gridGeometry_, curSol);
    }

    //! update all variables after grid adaption
    template<class SolutionVector>
    void updateAfterGridAdaption(const SolutionVector& curSol)
    {
        // update (always force data cache update as the grid changed)
        curGridVars_.init(*gridGeometry_, curSol);

        // for instationary problems also update the variables
        // for the previous time step to the new grid
        prevGridVars_ = curGridVars_;
    }

    /*!
     * \brief Sets the current state as the previous for next time step
     * \note this has to be called at the end of each time step
     */
    void advanceTimeStep()
    {
        prevGridVars_ = curGridVars_;
    }

    //! resets state to the one before time integration
    template<class SolutionVector>
    void resetTimeStep(const SolutionVector& solution)
    {
        // set the new time step variables to old variables
        curGridVars_ = prevGridVars_;
    }

    //! return the current variables
    const GridVariablesCache& curGridVars() const
    { return curGridVars_; }

    //! return the current variables
    GridVariablesCache& curGridVars()
    { return curGridVars_; }

    //! return the variables of the previous time step (for instationary problems)
    const GridVariablesCache& prevGridVars() const
    { return prevGridVars_; }

    //! return the variables of the previous time step (for instationary problems)
    GridVariablesCache& prevGridVars()
    { return prevGridVars_; }

    //! return the finite volume grid geometry
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

protected:
    std::shared_ptr<const GridGeometry> gridGeometry_; //!< pointer to the constant grid geometry

private:
    GridVariablesCache curGridVars_; //!< the current variables (primary and secondary variables)
    GridVariablesCache prevGridVars_; //!< the previous time step's variables (primary and secondary variables)
};

} // end namespace Dumux

#endif
