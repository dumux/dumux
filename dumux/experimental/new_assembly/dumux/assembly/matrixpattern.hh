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
 * \ingroup Assembly
 * \brief Classes and functions related to defining/modifying the
 *        non-zero entries in a sparse matrix pattern.
 */
#ifndef DUMUX_ASSEMBLY_MATRIX_PATTERN_HH
#define DUMUX_ASSEMBLY_MATRIX_PATTERN_HH

#include <algorithm>
#include <numeric>
#include <cassert>
#include <vector>
#include <set>

#include <dune/common/iteratorrange.hh>
#include <dune/common/iteratorfacades.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief Class to represent the non-zero entries in a matrix pattern.
 */
class MatrixPattern
{
    using RowIndices = std::set<std::size_t>;
    using IndexStorage = std::vector<RowIndices>;

    struct Index
    {
        std::size_t rowIndex;
        std::size_t colIndex;
    };

    class IndexIterator
    : public Dune::ForwardIteratorFacade<IndexIterator, const Index>
    {
    public:
        IndexIterator(const MatrixPattern& pattern, bool isEnd = false)
        : pattern_(pattern)
        , rowIndex_(isEnd ? pattern.indices_.size() : 0)
        {
            if (!isEnd_())
            {
                setBeginColIterator_();
                makeIndex_();
            }
        }

        void increment()
        {
            if (++colIt_; isEndColIterator_())
                if (++rowIndex_; !isEnd_())
                    setBeginColIterator_();

            if (!isEnd_())
                makeIndex_();
        }

        const Index& dereference() const
        { return index_; }

        bool equals(const IndexIterator& other) const
        {
            if (&pattern_ != &other.pattern_)
                return false;
            if (isEnd_() == other.isEnd_())
                return true;
            if (rowIndex_ != other.rowIndex_)
                return false;
            return colIt_ == other.colIt_;
        }

    private:
        bool isEnd_() const
        { return rowIndex_ == pattern_.indices_.size(); }

        bool isEndColIterator_() const
        { return colIt_ == pattern_.indices_[rowIndex_].end(); }

        void setBeginColIterator_()
        { colIt_ = pattern_.indices_[rowIndex_].begin(); }

        void makeIndex_()
        { index_ = Index{rowIndex_, *colIt_}; }

        const MatrixPattern& pattern_;
        std::size_t rowIndex_;
        typename RowIndices::iterator colIt_;
        Index index_;
    };

    friend class IndexIterator;

public:
    MatrixPattern() { resize(0, 0); }
    MatrixPattern(std::size_t rows, std::size_t cols) { resize(rows, cols); }

    void resize(std::size_t rows, std::size_t cols)
    {
        rows_ = rows;
        cols_ = cols;
        indices_.resize(rows);
    }

    void add(std::size_t i, std::size_t j)
    {
        assert(i < rows_);
        assert(j < cols_);
        indices_[i].insert(j);
    }

    auto begin() const { return IndexIterator{*this}; }
    auto end() const { return IndexIterator{*this, true}; }

    std::size_t numRows() const { return rows_; }
    std::size_t numCols() const { return cols_; }

    std::size_t nnz() const
    {
        return std::accumulate(
            indices_.begin(),
            indices_.end(),
            std::size_t{0},
            [] (std::size_t current, const auto& row) {
                return current + std::ranges::size(row);
            }
        );
    }

    std::size_t rowSize(std::size_t i) const
    {
        assert(i < indices_.size());
        return indices_[i].size();
    }

    friend auto indices(const MatrixPattern& pattern)
    {
        return Dune::IteratorRange<IndexIterator>(pattern.begin(),
                                                  pattern.end());
    }

private:
    std::size_t rows_;
    std::size_t cols_;
    std::vector<RowIndices> indices_;
};


/*!
 * \ingroup Assembly
 * \brief Turn the given matrix pattern into a flat pattern.
 * \tparam rowsPerBlock number of rows per entry of the given pattern.
 * \tparam colsPerBlock number of columns per entry of the given pattern.
 */
template<int rowsPerBlock, int colsPerBlock = rowsPerBlock>
MatrixPattern asFlatPattern(const MatrixPattern& pattern)
{
    MatrixPattern result(
        pattern.numRows()*rowsPerBlock,
        pattern.numCols()*colsPerBlock
    );

    for (const auto& index : indices(pattern))
    {
        const auto rowOffset = index.rowIndex*rowsPerBlock;
        const auto colOffset = index.colIndex*colsPerBlock;
        for (int i = 0; i < rowsPerBlock; ++i)
            for (int j = 0; j < colsPerBlock; ++j)
                result.add(rowOffset+i, colOffset+j);
    }

    return result;
}

} // namespace Dumux

#endif
