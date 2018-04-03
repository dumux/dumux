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
 * \brief The local element solution class for staggered methods
 */
#ifndef DUMUX_STAGGERED_ELEMENT_SOLUTION_HH
#define DUMUX_STAGGERED_ELEMENT_SOLUTION_HH

#include <type_traits>
#include <dune/istl/bvector.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux {

template<class PrimaryVariables>
using StaggeredElementSolution = Dune::BlockVector<PrimaryVariables>;

/*!
 * \ingroup StaggeredDiscretization
 * \brief  Make an element solution for staggered schemes
 * \note This is e.g. used to contruct an element solution at Dirichlet boundaries
 */
template<class FVElementGeometry, class PrimaryVariables>
auto elementSolution(PrimaryVariables&& priVars)
-> std::enable_if_t<FVElementGeometry::FVGridGeometry::discMethod == DiscretizationMethod::staggered,
                    StaggeredElementSolution<PrimaryVariables>>
{
    return StaggeredElementSolution<PrimaryVariables>({std::move(priVars)});
}

} // end namespace Dumux

#endif
