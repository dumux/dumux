// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Nonlinear
 * \brief A compressed subdomain-local view of a global Jacobian
 */
#ifndef DUMUX_NONLINEAR_PRECONDITIONING_MATRIX_VIEW_HH
#define DUMUX_NONLINEAR_PRECONDITIONING_MATRIX_VIEW_HH

#include <cstddef>
#include <vector>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>

#include <dumux/nonlinear/preconditioning/partition.hh>

namespace Dumux::NonlinearPreconditioning {

/*!
 * \ingroup Nonlinear
 * \brief A compressed, subdomain-owned view that the unmodified local assemblers can write into
 *
 * Rows are the degrees of freedom of subdomain \f$ N_i \f$; columns are either \f$ N_i \f$ (giving the
 * square block \f$ A_i = R_i J R_i^T \f$) or \f$ N_i \cup \Gamma_i \f$ (giving the rectangular block with
 * global columns that the exact RASPEN Jacobian requires).
 *
 * Writes addressed to a row outside \f$ N_i \f$ are discarded. This is the routine, expected path: an
 * assembly loop over \f$ N_i \cup \Gamma_i \f$ visits elements whose own neighbours lie further out, and
 * those rows are simply not wanted. Writes addressed to a column outside the column set are also
 * discarded, but that path should never be taken for a correctly constructed view, since every write's
 * column is the dof of the currently visited element.
 */
template<class MB>
class SubdomainMatrixView
{
public:
    using MatrixBlock = MB;
    using LocalMatrix = Dune::BCRSMatrix<MatrixBlock>;
    using Index = Partition::Index;
    using field_type = typename MatrixBlock::field_type;

    class RowProxy
    {
    public:
        RowProxy(SubdomainMatrixView& view, Index localRow)
        : view_(view), localRow_(localRow) {}

        MatrixBlock& operator[] (Index globalJ)
        {
            if (localRow_ == Partition::noIndex)
                return view_.sink_;

            const auto localCol = view_.colOf_[globalJ];
            if (localCol == Partition::noIndex)
            {
                ++view_.numColumnSinkWrites_;
                return view_.sink_;
            }

            return view_.matrix_[localRow_][localCol];
        }

    private:
        SubdomainMatrixView& view_;
        Index localRow_;
    };

    /*!
     * \brief Construct the view for one subdomain
     * \param partition the decomposition
     * \param influence the influence map the decomposition was built from
     * \param i the subdomain
     * \param includeFringeColumns whether to carry the \f$ \Gamma_i \f$ columns
     */
    SubdomainMatrixView(const Partition& partition,
                        const std::vector<std::vector<Index>>& influence,
                        SubdomainIndex i,
                        bool includeFringeColumns)
    : rowOf_(partition.numDofs(), Partition::noIndex)
    , colOf_(partition.numDofs(), Partition::noIndex)
    , rows_(partition.dofs(i))
    , columns_(includeFringeColumns ? partition.dofsWithFringe(i) : partition.dofs(i))
    , numRows_(partition.dofs(i).size())
    {
        for (Index k = 0; k < rows_.size(); ++k)
            rowOf_[rows_[k]] = k;

        const auto& columns = columns_;
        for (Index k = 0; k < columns.size(); ++k)
            colOf_[columns[k]] = k;
        numCols_ = columns.size();

        Dune::MatrixIndexSet pattern(numRows_, numCols_);
        for (const auto col : columns)
        {
            const auto localCol = colOf_[col];
            if (rowOf_[col] != Partition::noIndex)
                pattern.add(rowOf_[col], localCol);

            for (const auto row : influence[col])
                if (rowOf_[row] != Partition::noIndex)
                    pattern.add(rowOf_[row], localCol);
        }
        pattern.exportIdx(matrix_);
        matrix_ = 0.0;
    }

    RowProxy operator[] (Index globalI) { return RowProxy(*this, rowOf_[globalI]); }

    std::size_t N() const { return numRows_; }
    std::size_t M() const { return numCols_; }

    void setToZero()
    {
        matrix_ = 0.0;
        numColumnSinkWrites_ = 0;
    }

    const LocalMatrix& matrix() const { return matrix_; }
    LocalMatrix& matrix() { return matrix_; }

    //! the global dof of each local column, ascending
    const std::vector<Index>& columns() const { return columns_; }

    //! the global dof of each local row, ascending
    const std::vector<Index>& rows() const { return rows_; }

    //! must stay zero for a correctly constructed view; a non-zero count indicates a wrong loop range
    std::size_t numColumnSinkWrites() const { return numColumnSinkWrites_; }

private:
    LocalMatrix matrix_;
    std::vector<Index> rowOf_;
    std::vector<Index> colOf_;
    std::vector<Index> rows_;
    std::vector<Index> columns_;
    std::size_t numRows_;
    std::size_t numCols_;

    MatrixBlock sink_;
    std::size_t numColumnSinkWrites_ = 0;
};

} // end namespace Dumux::NonlinearPreconditioning

#endif
