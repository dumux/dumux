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
 * \brief MultiIndex access for dune block vector/matrix types.
 */
#ifndef DUMUX_LINEAR_DUNE_ACCESS_HH
#define DUMUX_LINEAR_DUNE_ACCESS_HH

#include <dune/common/indices.hh>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/common/multiindex.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {

template<typename V, typename Index> requires(
    Concepts::Indexable<V, Index>)
decltype(auto) getDuneVectorEntry(V& v, const Index& index)
{ return v[index]; }

template<typename V, Concepts::StaticIndex Index> requires(
    Concepts::Indexable<V, Dune::index_constant<Index::get()>>)
decltype(auto) getDuneVectorEntry(V& v, const Index& index)
{ return v[Dune::index_constant<Index::get()>()]; }

} // namespace Detail
#endif // DOXYGEN


struct DuneBlockVectorAccess
{
public:
    template<typename V, Concepts::MultiIndex I>
    static decltype(auto) get(V& v, const I& index)
    { return get_<0>(v, index); }

private:
    template<int i, typename V, Concepts::MultiIndex I>
    static decltype(auto) get_(V& v, const I& index)
    {
        if constexpr (i == I::size() - 1)
            return Detail::getDuneVectorEntry(v, index.template get<i>());
        else
            return get_<i+1>(
                Detail::getDuneVectorEntry(v, index.template get<i>()),
                index
            );
    }
};


struct DuneBlockMatrixAccess
{
public:
    template<typename M, Concepts::MultiIndex Row, Concepts::MultiIndex Col>
    static decltype(auto) get(M& m, const Row& row, const Col& col)
    { return get_<0, 0>(m, row, col); }

private:
    template<int i, int j, typename M, Concepts::MultiIndex Row, Concepts::MultiIndex Col>
    static decltype(auto) get_(M& m, const Row& row, const Col& col)
    {
        if constexpr (i == Row::size() - 1 && j == Col::size() - 1)
            return getEntry_(m, row.template get<i>(), col.template get<j>());
        else if constexpr (i == Row::size())
            return DuneBlockVectorAccess::get(m, getSubMultiIndex<j>(col));
        else if constexpr (j == Col::size())
            return DuneBlockVectorAccess::get(m, getSubMultiIndex<i>(row));
        else
            return get_<i+1, j+1>(
                getEntry_(m, row.template get<i>(), col.template get<j>()),
                row, col
            );
    }

    template<typename M, typename RowIndex, typename ColIndex>
    static decltype(auto) getEntry_(M& m,
                                    const RowIndex& row,
                                    const ColIndex& col)
    {
        return Detail::getDuneVectorEntry(
            Detail::getDuneVectorEntry(m, row),
            col
        );
    }
};

} // namespace Dumux

#endif
