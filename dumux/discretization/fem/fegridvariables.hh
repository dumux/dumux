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
 * \brief The grid variables class for finite element schemes. TODO Doc further
 */
#ifndef DUMUX_FE_GRID_VARIABLES_HH
#define DUMUX_FE_GRID_VARIABLES_HH

#include <type_traits>
#include <memory>
#include <cassert>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief The grid variables class for finite element schemes. TODO Doc further
 * \tparam the type of the grid geometry
 * \tparam the type of the secondary variables object
 */
template<class GG, class SV>
class FEGridVariables
{
public:
    //! export type of the finite volume grid geometry
    using GridGeometry = GG;

    //! export type of the volume variables
    using SecondaryVariables = SV;

    //! export primary variable type
    using PrimaryVariables = typename SecondaryVariables::PrimaryVariables;

    //! export scalar type (TODO get it directly from the volvars)
    using Scalar = std::decay_t<decltype(std::declval<PrimaryVariables>()[0])>;

    //! Constructor
    template<class Problem>
    FEGridVariables(std::shared_ptr<Problem> problem,
                    std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_(gridGeometry)
    {}

    //! initialize all variables (stationary case)
    template<class SolutionVector>
    void init(const SolutionVector& curSol)
    {}

    //! update all variables
    template<class SolutionVector>
    void update(const SolutionVector& curSol, bool forceCacheUpdate = false)
    {}

    //! update all variables after grid adaption
    template<class SolutionVector>
    void updateAfterGridAdaption(const SolutionVector& curSol)
    {}

    /*!
     * \brief Sets the current state as the previous for next time step
     */
    void advanceTimeStep()
    {}

    //! resets state to the one before time integration
    template<class SolutionVector>
    void resetTimeStep(const SolutionVector& solution)
    {}

    //! return the finite volume grid geometry
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

protected:
    std::shared_ptr<const GridGeometry> gridGeometry_; //!< pointer to the constant grid geometry
};

} // end namespace Dumux

#endif
