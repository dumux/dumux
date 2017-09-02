// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief A 2p1cni specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_MATRIX_CONVERTER
#define DUMUX_MATRIX_CONVERTER

#include "properties.hh"

#include <dumux/nonlinear/newtoncontroller.hh>
#include <dumux/linear/linearsolveracceptsmultitypematrix.hh>
#include "newtonconvergencewriter.hh"
#include <dune/common/tupleutility.hh>

namespace Dumux {

/*!
 * \ingroup PNMModel
 * \brief A PNM specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */

template <class MultiTypeBlockMatrix, class Scalar=double>
class MatrixConverter
{
    using MatrixBlock = typename Dune::FieldMatrix<Scalar, 1, 1>;
    using BCRSMatrix = typename Dune::BCRSMatrix<MatrixBlock>;

public:

    // copy the matrix and the vector to types the IterativeSolverBackend can handle
    static auto multiTypeToBCRSMatrix(const MultiTypeBlockMatrix &A)
    {
        // get the number of rows for the converted matrix
        const auto numRows = getNumRows_(A);

        auto M = BCRSMatrix(numRows, numRows, BCRSMatrix::random);

        setOccupationPattern_(M, A);
        copyValues_(M, A);

        return M;
    }

private:

    /*!
     * \brief Sets the occupation pattern (i.e. indices) for the converted matrix
     *
     * \param M The converted matrix
     * \param A The original subMatrix of the multitype blockmatrix
     */
    static void setOccupationPattern_(BCRSMatrix& M, const MultiTypeBlockMatrix& A)
    {
        // set the rowsizes
        const auto numRows = M.N();
        Dune::MatrixIndexSet occupationPattern;
        occupationPattern.resize(numRows, numRows);

        // lambda function to fill the occupation pattern
        auto addIndices = [&occupationPattern](const auto& subMatrix, const std::size_t startRow, const std::size_t startCol)
        {
            using BlockType = typename std::decay_t<decltype(subMatrix)>::block_type;
            const auto blockSizeI = BlockType::rows;
            const auto blockSizeJ = BlockType::cols;
            for(auto row = subMatrix.begin(); row != subMatrix.end(); ++row)
                for(auto col = row->begin(); col != row->end(); ++col)
                    for(std::size_t i = 0; i < blockSizeI; ++i)
                        for(std::size_t j = 0; j < blockSizeJ; ++j)
                            occupationPattern.add(startRow + row.index()*blockSizeI + i, startCol + col.index()*blockSizeJ + j);
        };

        int rowIndex = 0;
        Dune::Hybrid::forEach(A, [&addIndices, &rowIndex, &occupationPattern, numRows](const auto& rowOfMultiTypeMatrix)
        {
            int colIndex = 0;
            Dune::Hybrid::forEach(rowOfMultiTypeMatrix, [&addIndices, &occupationPattern, &colIndex, &rowIndex, numRows](const auto& subMatrix)
            {
                addIndices(subMatrix, rowIndex, colIndex);

                const auto numEq = std::decay_t<decltype(subMatrix)>::block_type::cols;
                colIndex += numEq*subMatrix.M();

                // if we arrived at the right side of the matrix, increase the row index
                if(colIndex == numRows)
                    rowIndex += numEq*subMatrix.N();
            });
        });

        occupationPattern.exportIdx(M);
    }

    /*!
     * \brief Sets the occupation pattern (i.e. indices) for the converted matrix
     *
     * \param M The converted matrix
     * \param A The original subMatrix of the multitype blockmatrix
     */
    static void copyValues_(BCRSMatrix& M, const MultiTypeBlockMatrix& A)
    {
        // set the rowsizes
        const auto numRows = M.N();

        // lambda function to copy the values
        auto copyValues = [&M](const auto& subMatrix, const std::size_t startRow, const std::size_t startCol)
        {
            using BlockType = typename std::decay_t<decltype(subMatrix)>::block_type;
            const auto blockSizeI = BlockType::rows;
            const auto blockSizeJ = BlockType::cols;
            for (auto row = subMatrix.begin(); row != subMatrix.end(); ++row)
                for (auto col = row->begin(); col != row->end(); ++col)
                    for (std::size_t i = 0; i < blockSizeI; ++i)
                        for (std::size_t j = 0; j < blockSizeJ; ++j)
                            M[startRow + row.index()*blockSizeI + i][startCol + col.index()*blockSizeJ + j] = subMatrix[row.index()][col.index()][i][j];

        };

        int rowIndex = 0;
        Dune::Hybrid::forEach(A, [&copyValues, &rowIndex, numRows](const auto& rowOfMultiTypeMatrix)
        {
            int colIndex = 0;
            Dune::Hybrid::forEach(rowOfMultiTypeMatrix, [&copyValues, &colIndex, &rowIndex, numRows](const auto& subMatrix)
            {
                copyValues(subMatrix, rowIndex, colIndex);

                const auto numEq = std::decay_t<decltype(subMatrix)>::block_type::cols;
                colIndex += numEq*subMatrix.M();

                // if we arrived at the right side of the matrix, increase the row index
                if(colIndex == numRows)
                    rowIndex += numEq*subMatrix.N();
            });
        });
    }

    /*!
     * \brief Calculates the total number of rows (== number of cols) for the converted matrix, assuming a block size of one
     *
     * \param A The original multitype blockmatrix
     */
    static std::size_t getNumRows_(const MultiTypeBlockMatrix& A)
    {
        // problem-agnostic version
        std::size_t numRows = 0;
        Dune::Hybrid::forEach(A[Dune::Indices::_0], [&numRows](const auto& subMatrixInFirstRow)
        {
            // std::cout << "in matrix of first row:\n";
            // std::cout << "numRows: " << subMatrixInFirstRow.N() << ", numCols: " << subMatrixInFirstRow.M() << std::endl;

            // The number of rows of the individual submatrice's block equals the respective number of equations.
            const auto numEq = std::decay_t<decltype(subMatrixInFirstRow)>::block_type::cols; // TODO: is this allways correct?
            numRows += numEq * subMatrixInFirstRow.M();
            // std::cout << "numEq: " << numEq << std::endl;
            // std::cout << "totalNumRows: " << numRows  << std::endl;
        });

        return numRows;
    }

};
}

#endif
