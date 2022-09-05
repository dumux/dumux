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
 * \brief Tests for accessing dune matrices by multi indices
 */
#include <iostream>

#include <dune/common/fmatrix.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>

#include <dumux/experimental/new_assembly/dumux/common/multiindex.hh>
#include <dumux/experimental/new_assembly/dumux/linear/dune/access.hh>


template<typename Accessor,
         typename BlockMatrix,
         typename RowIndexFactory,
         typename ColIndexFactory,
         typename DuneAccessor>
bool checkBlockMatrixAccess(const Accessor& acc,
                            const BlockMatrix& m,
                            int numRowBlocks,
                            int numColBlocks,
                            int blockRows,
                            int blockCols,
                            const RowIndexFactory& makeRowIndex,
                            const ColIndexFactory& makeColindex,
                            const DuneAccessor& duneAccessor)
{
    for (std::size_t i = 0; i < numRowBlocks; ++i)
        for (std::size_t j = 0; j < numColBlocks; ++j)
            for (int row = 0; row < blockRows; ++row)
                for (int col = 0; col < blockCols; ++col)
                    if (acc.get(m, makeRowIndex(i, row), makeColindex(j, col)) != duneAccessor(m, i, j, row, col))
                        return false;
    return true;
}

template<typename Accessor, typename BlockMatrix>
bool checkBlockMatrixAccess(const Accessor& acc,
                            const BlockMatrix& m,
                            int blockRows,
                            int blockCols)
{
    return checkBlockMatrixAccess(
        acc, m, m.N(), m.M(), blockRows, blockCols,
        [] (auto i0, auto i1) { return Dumux::MultiIndex{i0, i1}; },
        [] (auto i0, auto i1) { return Dumux::MultiIndex{i0, i1}; },
        [] (const auto& m, auto blockRow, auto blockCol, auto subRow, auto subCol) {
            return m[blockRow][blockCol][subRow][subCol];
        }
    );
}

template<int row, int col, typename Accessor, typename MultiTypeMatrix>
bool checkMultiTypeMatrixAccess(const Accessor& acc,
                                const MultiTypeMatrix& m,
                                int blockRows,
                                int blockCols)
{
    const auto& subMatrix = m[Dune::index_constant<row>()][Dune::index_constant<col>()];
    return checkBlockMatrixAccess(
        acc, m, subMatrix.N(), subMatrix.M(), blockRows, blockCols,
        [] (auto i0, auto i1) { return Dumux::MultiIndex{Dumux::_i<row>, i0, i1}; },
        [] (auto i0, auto i1) { return Dumux::MultiIndex{Dumux::_i<col>, i0, i1}; },
        [] (const auto& m, auto blockRow, auto blockCol, auto subRow, auto subCol) {
            return m[Dune::index_constant<row>()][Dune::index_constant<col>()][blockRow][blockCol][subRow][subCol];
        }
    );
}

int main (int argc, char *argv[])
{
    using namespace Dumux;

    DuneBlockMatrixAccess accessor;
    {  // dynamic access - 2 levels
        Dune::DynamicMatrix<Dune::FieldMatrix<int, 2, 3>> m{{
            {
                {{41, 42, 43}, {44, 45, 46}},
                {{47, 48, 49}, {50, 51, 52}}
            },
            {
                {{53, 54, 55}, {56, 57, 58}},
                {{59, 60, 61}, {62, 63, 64}}
            }
        }};

        if (!checkBlockMatrixAccess(accessor, m, 2, 3))
            DUNE_THROW(Dune::InvalidStateException, "Unexpected block");
    }

    {  // static access
        constexpr int size0 = 2;
        constexpr int size1 = 1;
        using BT00 = Dune::DynamicMatrix<Dune::FieldMatrix<int, size0, size0>>;
        using BT01 = Dune::DynamicMatrix<Dune::FieldMatrix<int, size0, size1>>;
        using BT10 = Dune::DynamicMatrix<Dune::FieldMatrix<int, size1, size0>>;
        using BT11 = Dune::DynamicMatrix<Dune::FieldMatrix<int, size1, size1>>;
        using Row0 = Dune::MultiTypeBlockVector<BT00, BT01>;
        using Row1 = Dune::MultiTypeBlockVector<BT10, BT11>;

        Row0 row0;
        row0[Dune::Indices::_0] = BT00{{
            {  // row0
                {  // matrix00
                    {{41, 42}},
                    {{43, 44}}
                },
                {  // matrix01
                    {{45, 46}},
                    {{47, 48}}
                }
            },
            {  // row1
                {  // matrix10
                    {{49, 50}},
                    {{51, 52}}
                },
                {  // matrix11
                    {{53, 54}},
                    {{55, 56}}
                }
            }
        }};
        row0[Dune::Indices::_1] = BT01{{
            {  // row0
                {  // matrix00
                    {{57}},
                    {{58}},
                },
                {  // matrix01
                    {{59}},
                    {{60}}
                }
            },
            {  // row1
                {  // matrix10
                    {{61}},
                    {{62}}
                },
                {  // matrix11
                    {{63}},
                    {{64}}
                }
            }
        }};

        Row1 row1;
        row1[Dune::Indices::_0] = BT10{{
            {  // row0
                {  // matrix00
                    {{65, 66}}
                },
                {  // matrix01
                    {{67, 68}}
                }
            },
            {  // row1
                {  // matrix10
                    {{69}, {70}}
                },
                {  // matrix11
                    {{71}, {72}}
                }
            }
        }};
        row1[Dune::Indices::_1] = BT11{{
            {  // row0
                {  // matrix00
                    {{{73}}}
                },
                {  // matrix01
                    {{{74}}}
                }
            },
            {  // row1
                {  // matrix10
                    {{{75}}}
                },
                {  // matrix11
                    {{{76}}}
                }
            }
        }};

        Dune::MultiTypeBlockMatrix<Row0, Row1> m;
        m[Dune::Indices::_0] = row0;
        m[Dune::Indices::_1] = row1;

        DuneBlockMatrixAccess accessor;

        // upper left block
        if (!checkMultiTypeMatrixAccess<0, 0>(accessor, m, size0, size0))
            DUNE_THROW(Dune::InvalidStateException, "Unexpected block");

        // upper right block
        if (!checkMultiTypeMatrixAccess<0, 1>(accessor, m, size0, size1))
            DUNE_THROW(Dune::InvalidStateException, "Unexpected block");

        // lower left block
        if (!checkMultiTypeMatrixAccess<1, 0>(accessor, m, size1, size0))
            DUNE_THROW(Dune::InvalidStateException, "Unexpected block");

        // lower right block
        if (!checkMultiTypeMatrixAccess<1, 1>(accessor, m, size1, size1))
            DUNE_THROW(Dune::InvalidStateException, "Unexpected block");
    }

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
