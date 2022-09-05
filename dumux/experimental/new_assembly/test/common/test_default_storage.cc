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
 * \brief Tests for the size helpers
 */
#include <iostream>
#include <type_traits>

#include <dune/common/exceptions.hh>

#include <dumux/experimental/new_assembly/dumux/common/typetraits.hh>
#include <dumux/experimental/new_assembly/dumux/common/storage.hh>

int main (int argc, char *argv[])
{
    Dumux::DefaultStorage<int, 3> staticFromInt{};
    Dumux::DefaultStorage<int, Dumux::Size::dynamic> fromDynamic{};

    static_assert(std::is_same_v<decltype(staticFromInt), Dune::ReservedVector<int, 3>>);
    static_assert(std::is_same_v<decltype(fromDynamic), std::vector<int>>);

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
