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
 * \ingroup Common
 * \brief Defines the index types used for grid and local indices.
 */
#ifndef DUMUX_COMMON_INDEX_TRAITS_HH
#define DUMUX_COMMON_INDEX_TRAITS_HH

#include <cstdint>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief Struture to define the index types used for grid and local indices.
 * \tparam GridView The grid view type
 */
template<class GridView>
struct IndexTraits
{
    using GridIndex = typename GridView::IndexSet::IndexType;
    using LocalIndex = unsigned int;
    using SmallLocalIndex = std::uint_least8_t;
};

} // namespace Dumux

#endif
