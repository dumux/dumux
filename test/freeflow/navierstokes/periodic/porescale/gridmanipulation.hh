// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Grid Manipulation
 * \brief Manipulates host, sub, extended and convoluted grids
 */

#ifndef DUMUX_GRID_MANIPULATION_HH
#define DUMUX_GRID_MANIPULATION_HH

#include <tuple>
#include <cmath>
#include <limits>
#include <sstream>

#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/entitymap.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectingentities.hh>

namespace Dumux::GridManipulation {

using Scalar = double;
static const int dimWorld = 2;
using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
using VectorQuantity = GlobalPosition;
using IntArray = std::array<int, dimWorld>;

using Geometry = Dune::MultiLinearGeometry<Scalar, dimWorld, dimWorld>;
using GeoEntitySet = GeometriesEntitySet<Geometry>;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// Tile and Trim for Conv Filter /////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class Image, class SubGridView>
std::vector<int> colsPerRowInImage(const Image& img,
                                   const SubGridView& subGridView)
{
    std::vector<int> colsPerRow(img.header().nRows, 0);
    Scalar rowHeight = std::numeric_limits<double>::lowest();
    int rowIdx = -1;
    for (const auto& element : elements(subGridView))
    {
        Scalar height = element.geometry().center()[1];
        if (height > rowHeight)
        {
            rowHeight = height;
            rowIdx++;
        }
        colsPerRow[rowIdx]++;
    }
    return colsPerRow;
}

template<class Content>
std::vector<std::vector<Content>> wrapVectorToMatrix(const std::vector<Content>& originalVector,
                                                     const IntArray& matrixDims)
{
    int numCols = matrixDims[0];
    int numRows = matrixDims[1];
    std::vector<std::vector<Content>> wrappedMatrix(numRows, std::vector<Content>(numCols));

    int k = 0;
    for (int i = 0; i < numRows; i++)
    {
        for (int j = 0; j < numCols; j++)
        {
            wrappedMatrix[i][j] = originalVector[k];
            k++;
        }
    }

    return wrappedMatrix;
}

template <class Content>
std::vector<Content> unwrapSquareMatrix(const std::vector<std::vector<Content>>& originalMatrix)
{
    std::vector<Content> unwrappedVector(originalMatrix.size() * originalMatrix[0].size());
    int k = 0;
    for (int i = 0; i < originalMatrix.size() ; i++)
    {
        for (int j = 0; j < originalMatrix[0].size(); j++)
        {
            unwrappedVector[k] = originalMatrix[i][j];
            k++;
        }
    }

    return unwrappedVector;
}

std::string positionToString(const GlobalPosition& pos, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    for (const auto& coord : pos)
        out << std::fixed << coord;
    return out.str();
}

// Concatenates two matrices A and B horizontally to create a new matrix C
template<class Content>
std::vector<std::vector<Content>> concatenateMatrices(const std::vector<std::vector<Content>>& A,
                                                      const std::vector<std::vector<Content>>& B)
{
    std::vector<std::vector<Content>> outputMatrix;

    // Copy the elements of A and B into C
    for (int i = 0; i < A.size(); i++)
    {
        std::vector<Content> row;
        for (const auto& element : A[i])
            row.push_back(element);
        for (const auto& element : B[i])
            row.push_back(element);
        outputMatrix.push_back(row);
    }

    return outputMatrix;
}

// Stacks a matrix A vertically to create a new matrix
template<class Content>
std::vector<std::vector<Content>> stackMatrixVertically(const std::vector<std::vector<Content>>& matrix,
                                                        const int& x)
{
    std::vector<std::vector<Content>> outputMatrix;

    // Copy the elements of matrix into C x times
    for (int i = 0; i < x; i++)
        for (const auto& row : matrix)
            outputMatrix.push_back(row);

    return outputMatrix;
}

template<class Content>
std::vector<std::vector<Content>> tileVelocityMatrix(const std::vector<std::vector<Content>>& wrappedVelocityMatrix)
{
    // Develop the same matrix tiled to a 3x3
    std::vector<std::vector<Content>> concatenatedMatrix = concatenateMatrices(wrappedVelocityMatrix, concatenateMatrices(wrappedVelocityMatrix, wrappedVelocityMatrix));
    return stackMatrixVertically(concatenatedMatrix, 3);
}

template<class Content>
std::vector<std::vector<Content>> tileTracerMatrix(const std::vector<std::vector<Content>>& wrappedInputTracerMatrix,
                                                   const std::vector<std::vector<Content>>& wrappedTracerMatrix)
{
    // Develop the same matrix tiled to a 3x3
    std::vector<std::vector<Content>> concatenatedMatrix = concatenateMatrices(wrappedInputTracerMatrix, wrappedTracerMatrix);
    return stackMatrixVertically(concatenatedMatrix, 3);
}

template <class Content>
std::vector<Content> compileTracerSolution(const std::vector<Content>& tracerExtended)
{
    IntArray tripleGridDims = {600, 200};
    std::vector<std::vector<Content>> wrappedTracerMatrix = wrapVectorToMatrix(tracerExtended, tripleGridDims);
    std::vector<std::vector<Content>> tiledMatrix = stackMatrixVertically(wrappedTracerMatrix, 3);
    return GridManipulation::unwrapSquareMatrix(tiledMatrix);
}

template <class Content>
std::vector<Content> compileVelocitySolution(const std::vector<Content>& velocitySolExtended)
{
    IntArray singleGridDims = {200, 200};
    std::vector<std::vector<Content>> wrappedVelocityMatrix = wrapVectorToMatrix(velocitySolExtended, singleGridDims);
    std::vector<std::vector<Content>> tiledMatrix = tileVelocityMatrix(wrappedVelocityMatrix);
    return GridManipulation::unwrapSquareMatrix(tiledMatrix);
}


template<class Content>
std::vector<Content> trimCompiledSolution(const std::vector<Content>& originalVector,
                                          const std::array<int, 2>& filterSize)
{
    IntArray tiledGridDimensions = {600, 600};
    std::vector<std::vector<Content>> wrappedTiledTracerMatrix = wrapVectorToMatrix(originalVector, tiledGridDimensions);

    int rows = wrappedTiledTracerMatrix.size()/3;
    int cols = wrappedTiledTracerMatrix[0].size()/3;

    std::array<int, 2> cellExtension;
    cellExtension[0] = (filterSize[0] - 1) / 2;
    cellExtension[1] = (filterSize[1] - 1) / 2;

    std::array<int, 2> lowerLeftIndexes = {(rows - cellExtension[0]), (cols - cellExtension[1])};
    std::array<int, 2> upperRightIndexes = {(2 * rows + cellExtension[0]), (2 * cols + cellExtension[1])};

    std::vector<std::vector<Content>> trimmedMatrix(upperRightIndexes[0] - lowerLeftIndexes[0],
                                                    std::vector<Content>(upperRightIndexes[1] - lowerLeftIndexes[1]));

    for (int i = lowerLeftIndexes[0]; i < upperRightIndexes[0]; i++)
        for (int j = lowerLeftIndexes[1]; j < upperRightIndexes[1]; j++)
            trimmedMatrix[i - lowerLeftIndexes[0]][j - lowerLeftIndexes[1]] = wrappedTiledTracerMatrix[i][j];

    return GridManipulation::unwrapSquareMatrix(trimmedMatrix);
}

template<class Content>
std::vector<Content> trimCompiledSolution(const std::vector<Content>& originalVector,
                                          const Scalar& filterRadius)
{
    IntArray tiledGridDimensions = {600, 600};
    std::vector<std::vector<Content>> wrappedTiledTracerMatrix = wrapVectorToMatrix(originalVector, tiledGridDimensions);

    int rows = wrappedTiledTracerMatrix.size()/3;
    int cols = wrappedTiledTracerMatrix[0].size()/3;

    std::array<int, 2> cellExtension;
    cellExtension[0] = rows*filterRadius;
    cellExtension[1] = cols*filterRadius;

    std::array<int, 2> lowerLeftIndexes = {(rows - cellExtension[0]), (cols - cellExtension[1])};
    std::array<int, 2> upperRightIndexes = {(2 * rows + cellExtension[0]), (2 * cols + cellExtension[1])};

    std::vector<std::vector<Content>> trimmedMatrix(upperRightIndexes[0] - lowerLeftIndexes[0],
                                                    std::vector<Content>(upperRightIndexes[1] - lowerLeftIndexes[1]));

    for (int i = lowerLeftIndexes[0]; i < upperRightIndexes[0]; i++)
        for (int j = lowerLeftIndexes[1]; j < upperRightIndexes[1]; j++)
            trimmedMatrix[i - lowerLeftIndexes[0]][j - lowerLeftIndexes[1]] = wrappedTiledTracerMatrix[i][j];

    return GridManipulation::unwrapSquareMatrix(trimmedMatrix);
}

}  // end namespace Dumux::GridManipulation

#endif
