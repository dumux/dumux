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
#include <concepts>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/experimental/new_assembly/dumux/common/typetraits.hh>

struct Incomplete;
struct Complete {};
struct Indexable { int operator[](std::size_t i) { return 0; } };


int main (int argc, char *argv[])
{
    using namespace Dumux;

    static_assert(!isComplete<Incomplete>);
    static_assert(isIncomplete<Incomplete>);

    static_assert(isComplete<Complete>);
    static_assert(!isIncomplete<Complete>);

    static_assert(isIndexable<std::vector<Complete>>);
    static_assert(isIndexable<std::vector<Complete>, int>);

    static_assert(
        std::is_same_v<
            IndexedType<std::vector<Complete>>,
            Complete
        >
    );

    static_assert(indexDepth<std::vector<int>> == 1);
    static_assert(indexDepth<std::vector<std::vector<int>>> == 2);
    static_assert(indexDepth<std::vector<Indexable>> == 2);

    if (DynamicSize{}*DynamicSize{} != DynamicSize{})
        DUNE_THROW(Dune::InvalidStateException, "Unexpected operator*");
    if (DynamicSize{}+1 != DynamicSize{})
        DUNE_THROW(Dune::InvalidStateException, "Unexpected operator+");
    if (DynamicSize{}*1 != DynamicSize{})
        DUNE_THROW(Dune::InvalidStateException, "Unexpected operator*int{}");
    if (Size::pow<DynamicSize{}, 2>() != DynamicSize{})
        DUNE_THROW(Dune::InvalidStateException, "Unexpected power result");

    if (DynamicSize{} != DynamicSize{})
        DUNE_THROW(Dune::InvalidStateException, "Expected dynamic sizes to always compare equal");
    if (!(DynamicSize{} == DynamicSize{}))
        DUNE_THROW(Dune::InvalidStateException, "Expected dynamic sizes to always compare equal");
    if (DynamicSize{} < DynamicSize{})
        DUNE_THROW(Dune::InvalidStateException, "Unexpected operator<");
    if (DynamicSize{} > DynamicSize{})
        DUNE_THROW(Dune::InvalidStateException, "Unexpected operator>");

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
