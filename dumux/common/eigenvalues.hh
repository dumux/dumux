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
 * \brief Defines a functions to calculate eigenvalues and eigenvectors of n x n matrices. For n > 2 a cyclic jacobi method is used.
 * This implementation is not efficient for larg matrices!
 */
#ifndef DUMUX_EIGENVALUES_HH
#define DUMUX_EIGENVALUES_HH

#include <algorithm>
#include <cmath>

#include "math.hh"

namespace Dumux
{

template<int dim, class Matrix>
void identityMatrix(Matrix& matrix)
{
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            matrix[i][j] = (i == j) ? 1.0 : 0.0;
}

template<int dim, class Matrix>
double calcOffDiagonalNorm(Matrix& matrix)
{
    double norm = 0;
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            if (i != j)
            {
                norm += matrix[i][j] * matrix[i][j];
            }

    using std::sqrt;
    return sqrt(norm);
}

/*!
 * \briefFunction to calculate eigenvalues of n x n matrices
 *
 * \param eigVel Vector for storing the eigenvalues
 * \param matrix n x n matrices for which eigenvalues have to be calculated
 * \param relativeTolerance tolerance for the relative convergence criterion (default: 0.01)
 */
template<int dim, class EVVectorType, class MatrixType>
bool calculateEigenValues(EVVectorType &eigVel, MatrixType& matrix, double relativeTolerance = 0.01)
{
    if (dim == 2)
    {
        eigVel = 0;

        double b = -(matrix[0][0] + matrix[1][1]);
        double c = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

        eigVel[0] = (-b + sqrt(b * b - 4.0 * c)) / 2.0;
        eigVel[1] = (-b - sqrt(b * b - 4.0 * c)) / 2.0;

        using std::isnan;
        using std::isinf;
        if (isnan(eigVel[0]) || isinf(eigVel[0]))
            return false;

        if (isnan(eigVel[1]) || isinf(eigVel[1]))
            return false;

        return true;
    }
    else if (dim > 2)
    {
        int maxIter = 100;
        int iter = 0;
        double matrixNorm = matrix.frobenius_norm();
        double offDiagonalNorm = calcOffDiagonalNorm<dim>(matrix);
        MatrixType rotationMatrix(0.0);
        MatrixType evMatrix(matrix);

        while (iter < maxIter && offDiagonalNorm > relativeTolerance * matrixNorm)
        {
            for (int i = 0; i < dim - 1; i++)
            {
                for (int j = i + 1; j < dim; j++)
                {
                    identityMatrix<dim>(rotationMatrix);

                    double theta = (evMatrix[i][i] - evMatrix[j][j])
                            / (2 * evMatrix[i][j]);
                    using std::abs;
                    using std::sqrt;
                    double t = sign(theta) / (abs(theta) + sqrt(1 + theta * theta));
                    double c = 1 / sqrt(1 + t * t);
                    double s = c * t;

                    rotationMatrix[i][i] = c;
                    rotationMatrix[j][j] = c;
                    rotationMatrix[i][j] = s;
                    rotationMatrix[j][i] = -s;

                    evMatrix.leftmultiply(rotationMatrix);

                    rotationMatrix[i][j] = -s;
                    rotationMatrix[j][i] = s;

                    evMatrix.rightmultiply(rotationMatrix);
                }
            }
            matrixNorm = evMatrix.frobenius_norm();
            offDiagonalNorm = calcOffDiagonalNorm<dim>(evMatrix);
            iter++;
        }

        for (int i = 0; i < dim; i++)
        {
            eigVel[i] = evMatrix[i][i];
            using std::isinf;
            using std::isnan;
            if (isnan(eigVel[i]) || isinf(eigVel[i]))
                return false;
        }

        return true;
    }

    return false;
}

//! Function to calculate eigenvalues and eigenvectors of n x n matrices
/*
 * \param eigVel Vector for storing the eigenvalues
 * \param eigVec n x n for storing the eigenvectors
 * \param matrix n x n matrices for which eigenvalues have to be calculated
 * \param relativeTolerance tolerance for the relative convergence criterion (default: 0.01)
 */
template<int dim, class EVVectorType, class MatrixType>
bool calculateEigenValues(EVVectorType &eigVel, MatrixType& eigVec, MatrixType& matrix, double relativeTolerance = 0.01)
{
        int maxIter = 100;
        int iter = 0;
        double matrixNorm = matrix.frobenius_norm();
        double offDiagonalNorm = calcOffDiagonalNorm<dim>(matrix);
        MatrixType rotationMatrix(0.0);
        MatrixType evMatrix(matrix);
        MatrixType evecMatrix(matrix);
        identityMatrix<dim>(evecMatrix);

        while (iter < maxIter && offDiagonalNorm > relativeTolerance * matrixNorm)
        {
            for (int i = 0; i < dim - 1; i++)
            {
                for (int j = i + 1; j < dim; j++)
                {
                    identityMatrix<dim>(rotationMatrix);

                    double theta = (evMatrix[i][i] - evMatrix[j][j])
                            / (2 * evMatrix[i][j]);
                    using std::abs;
                    using std::sqrt;
                    double t = sign(theta) / (abs(theta) + sqrt(1 + theta * theta));
                    double c = 1 / sqrt(1 + t * t);
                    double s = c * t;

                    rotationMatrix[i][i] = c;
                    rotationMatrix[j][j] = c;
                    rotationMatrix[i][j] = s;
                    rotationMatrix[j][i] = -s;

                    evMatrix.leftmultiply(rotationMatrix);

                    rotationMatrix[i][j] = -s;
                    rotationMatrix[j][i] = s;

                    evMatrix.rightmultiply(rotationMatrix);
                    evecMatrix.rightmultiply(rotationMatrix);
                }
            }
            matrixNorm = evMatrix.frobenius_norm();
            offDiagonalNorm = calcOffDiagonalNorm<dim>(evMatrix);
            iter++;
        }

        for (int i = 0; i < dim; i++)
        {
            eigVel[i] = evMatrix[i][i];
            using std::isinf;
            using std::isnan;
            if (isnan(eigVel[i]) || isinf(eigVel[i]))
                return false;
            for (int j = 0; j < dim; j++)
            {
                eigVec[i][j] = evecMatrix[j][i];
            }
        }

        return true;
}


}

#endif
