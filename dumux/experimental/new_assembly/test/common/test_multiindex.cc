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
#include <algorithm>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/experimental/new_assembly/dumux/common/multiindex.hh>

int main (int argc, char *argv[])
{
    using namespace Dumux;
    using namespace Dumux::Indices;

    constexpr MultiIndex index{1, _0, 2};
    static_assert(index.template get<0>() == 1);
    static_assert(index.template get<1>().get() == 0);
    static_assert(index.template get<2>() == 2);
    std::cout << index << std::endl;

    constexpr MultiIndex index2{1, _1, index};
    static_assert(index2.template get<0>() == 1);
    static_assert(index2.template get<1>().get() == 1);
    static_assert(index2.template get<2>().template get<0>() == 1);
    static_assert(index2.template get<2>().template get<1>().get() == 0);
    static_assert(index2.template get<2>().template get<2>() == 2);
    std::cout << index2 << std::endl;

    constexpr auto subIndex = getSubMultiIndex<1, 2>(index2);
    static_assert(subIndex.template get<0>().get() == 1);
    static_assert(subIndex.template get<1>().template get<0>() == 1);
    static_assert(subIndex.template get<1>().template get<1>().get() == 0);
    static_assert(subIndex.template get<1>().template get<2>() == 2);
    std::cout << subIndex << std::endl;

    constexpr auto subIndex2 = getSubMultiIndex<0, 1>(index2);
    static_assert(subIndex2.template get<0>() == 1);
    static_assert(subIndex2.template get<1>().get() == 1);
    std::cout << subIndex2 << std::endl;

    if (!std::ranges::equal(asArray(index), std::vector<int>{{1, 0, 2}}))
        DUNE_THROW(Dune::InvalidStateException, "Conversion to array failed");

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
