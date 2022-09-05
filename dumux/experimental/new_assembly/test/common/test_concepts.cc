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
 * \brief Tests for the MultiIndex class
 */
#include <iostream>
#include <vector>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>

struct IndexableType { int operator[](std::size_t i) { return 0; } };

int main()
{
    using namespace Dumux::Concepts;

    static_assert(MDArray<IndexableType>);
    static_assert(MDArray<IndexableType, 1>);
    static_assert(MDArray<std::vector<IndexableType>, 2>);
    static_assert(!MDArray<std::vector<IndexableType>, 1>);
    static_assert(!MDArray<std::vector<IndexableType>, 3>);

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
