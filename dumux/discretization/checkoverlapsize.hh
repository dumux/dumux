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
 * \brief Check the overlap size for different discretization methods
 */
#ifndef DUMUX_DISCRETIZATION_CHECK_OVERLAP_SIZE_HH
#define DUMUX_DISCRETIZATION_CHECK_OVERLAP_SIZE_HH

#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Check if the overlap size is valid for a given discretization method
 * \note the default checks if the grid has at least an overlap of one if there are no ghosts
 * \note for sequential grids every overlap is fine
 * \note specialize this for your discretization method if the default doesn't apply
 */
template<DiscretizationMethod discMethod>
struct CheckOverlapSize
{
    template<class GridView>
    static bool isValid(const GridView& gridView) noexcept
    { return gridView.comm().size() <= 1 || gridView.overlapSize(0) + gridView.ghostSize(0) > 0; }
};

//! specialization for the box method which requires an overlap size of 0
template<>
struct CheckOverlapSize<DiscretizationMethod::box>
{
    template<class GridView>
    static bool isValid(const GridView& gridView) noexcept
    { return gridView.comm().size() <= 1 || gridView.overlapSize(0) == 0; }
};

//! specialization for the finite element method which requires an overlap size of 0
//! \note Overloads for bases that require overlap regions can be defined in the future
template<>
struct CheckOverlapSize<DiscretizationMethod::fem>
{
    template<class FEBasis>
    static bool isValid(const FEBasis& feBasis) noexcept
    { return feBasis.gridView().comm().size() <= 1 || feBasis.gridView().overlapSize(0) == 0; }
};

} // end namespace Dumux

#endif
