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

#include <dune/common/exceptions.hh>

#include <dumux/experimental/new_assembly/dumux/common/multiindex.hh>
#include <dumux/experimental/new_assembly/dumux/common/indexstrategies.hh>

int main (int argc, char *argv[])
{
    using namespace Dumux;
    using namespace Dumux::Indices;

    {  // blocked strategy
        BlockedIndexStrategy<2> blockedStrategy;
        MultiIndex index{1, _0};
        index = blockedStrategy[index];

        if (index.template get<0>() != 1)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected index");
        if (index.template get<1>().get() != 0)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected index");
    }

    {  // flat strategy (depth 2)
        FlatIndexStrategy<2> flatStrategy{10, 2};

        auto mapped = flatStrategy[MultiIndex{5, 0}];
        if (mapped.size() != 1)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected flat index depth");
        if (mapped.template get<0>() != 10)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected flat index");

        mapped = flatStrategy[MultiIndex{9, 1}];
        if (mapped.size() != 1)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected flat index depth");
        if (mapped.template get<0>() != 19)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected flat index");
    }

    {  // flat strategy (depth 3)
        FlatIndexStrategy<3> flatStrategy{_1, 10, 2};

        auto mapped = flatStrategy[MultiIndex{_0, 5, 0}];
        if (mapped.size() != 1)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected flat index depth");
        if (mapped.template get<0>() != 10)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected flat index");

        mapped = flatStrategy[MultiIndex{_1, 5, 0}];
        if (mapped.size() != 1)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected flat index depth");
        if (mapped.template get<0>() != 30)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected flat index");
    }

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
