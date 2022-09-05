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

#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dumux/experimental/new_assembly/dumux/common/multiindex.hh>
#include <dumux/experimental/new_assembly/dumux/linear/dune/access.hh>

int main (int argc, char *argv[])
{
    using namespace Dumux;

    DuneBlockVectorAccess accessor;
    {  // dynamic access - 2 levels
        Dune::DynamicVector<Dune::FieldVector<int, 2>> v{{
            {41, 42},
            {43, 44}
        }};
        if (accessor.get(v, MultiIndex{0, 0}) != 41)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
        if (accessor.get(v, MultiIndex{0, 1}) != 42)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
        if (accessor.get(v, MultiIndex{1, 0}) != 43)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
        if (accessor.get(v, MultiIndex{1, 1}) != 44)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
    }

    {  // dynamic access - 3 levels
        Dune::DynamicVector<Dune::FieldVector<Dune::FieldVector<int, 2>, 2>> v{{
            {{41, 42}, {43, 44}},
            {{45, 46}, {47, 48}}
        }};
        if (accessor.get(v, MultiIndex{0, 0, 0}) != 41)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
        if (accessor.get(v, MultiIndex{0, 0, 1}) != 42)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
        if (accessor.get(v, MultiIndex{0, 1, 0}) != 43)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
        if (accessor.get(v, MultiIndex{0, 1, 1}) != 44)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");

        if (accessor.get(v, MultiIndex{1, 0, 0}) != 45)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
        if (accessor.get(v, MultiIndex{1, 0, 1}) != 46)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
        if (accessor.get(v, MultiIndex{1, 1, 0}) != 47)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
        if (accessor.get(v, MultiIndex{1, 1, 1}) != 48)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
    }

    {  // static access
        Dune::MultiTypeBlockVector<
            Dune::DynamicVector<Dune::FieldVector<int, 1>>,
            Dune::DynamicVector<Dune::FieldVector<int, 2>>
        > v;
        v[Dune::Indices::_0] = Dune::DynamicVector<Dune::FieldVector<int, 1>>{{42}};
        v[Dune::Indices::_1] = Dune::DynamicVector<Dune::FieldVector<int, 2>>{{
            {43, 44},
            {45, 46}
        }};

        if (accessor.get(v, MultiIndex{_0, 0, 0}) != 42)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
        if (accessor.get(v, MultiIndex{_1, 0, 0}) != 43)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
        if (accessor.get(v, MultiIndex{_1, 0, 1}) != 44)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
        if (accessor.get(v, MultiIndex{_1, 1, 0}) != 45)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
        if (accessor.get(v, MultiIndex{_1, 1, 1}) != 46)
            DUNE_THROW(Dune::InvalidStateException, "Unexpected value");
    }

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
