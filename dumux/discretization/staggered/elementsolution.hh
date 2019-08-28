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
 * \ingroup StaggeredDiscretization
 * \brief The local element solution class for staggered methods
 */
#ifndef DUMUX_STAGGERED_ELEMENT_SOLUTION_HH
#define DUMUX_STAGGERED_ELEMENT_SOLUTION_HH

#include <type_traits>
#include <dune/istl/bvector.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief  Helper function to create a PrimaryVariables object from CellCenterPrimaryVariables
 * \tparam PrimaryVariables The type of the desired primary variables object
 * \tparam CellCenterPrimaryVariables The type of the cell center (input) primary variables object
 * \param cellCenterPriVars The cell center (input) primary variables object
 */
template<class PrimaryVariables, class CellCenterPrimaryVariables>
PrimaryVariables makePriVarsFromCellCenterPriVars(const CellCenterPrimaryVariables& cellCenterPriVars)
{
    static_assert(int(PrimaryVariables::dimension) > int(CellCenterPrimaryVariables::dimension),
                  "PrimaryVariables' size must be greater than the one of CellCenterPrimaryVariables");

    PrimaryVariables priVars(0.0);
    constexpr auto offset = PrimaryVariables::dimension - CellCenterPrimaryVariables::dimension;
    for (std::size_t i = 0; i < cellCenterPriVars.size(); ++i)
        priVars[i + offset] = cellCenterPriVars[i];
    return priVars;
}

template<class PrimaryVariables>
using StaggeredElementSolution = Dune::BlockVector<PrimaryVariables>;

/*!
 * \ingroup StaggeredDiscretization
 * \brief  Make an element solution for staggered schemes
 * \note This is e.g. used to construct an element solution at Dirichlet boundaries
 */
template<class FVElementGeometry, class PrimaryVariables>
auto elementSolution(PrimaryVariables&& priVars)
-> std::enable_if_t<FVElementGeometry::GridGeometry::discMethod == DiscretizationMethod::staggered,
                    StaggeredElementSolution<PrimaryVariables>>
{
    return StaggeredElementSolution<PrimaryVariables>({std::move(priVars)});
}

/*!
 * \ingroup StaggeredDiscretization
 * \brief  Helper function to create an elementSolution from cell center primary variables
 * \tparam PrimaryVariables The type of the desired primary variables object
 * \tparam CellCenterPrimaryVariables The type of the cell center (input) primary variables object
 * \param cellCenterPriVars The cell center (input) primary variables object
 */
template<class PrimaryVariables, class CellCenterPrimaryVariables>
StaggeredElementSolution<PrimaryVariables> makeElementSolutionFromCellCenterPrivars(const CellCenterPrimaryVariables& cellCenterPriVars)
{
    return StaggeredElementSolution<PrimaryVariables>({makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars)});
}

} // end namespace Dumux

#endif
