// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 * \brief Define some often used mathematical functions
 */
#ifndef DUNE_MATH_HH
#define DUNE_MATH_HH

namespace Dune
{

/*!
 * \brief Calculate the harmonic mean of two scalar values.
 */
template <class Scalar>
Scalar harmonicMean(Scalar x, Scalar y)
{
    if (x == 0 || y == 0)
        return 0;
    return (2*x*y)/(x + y);
};

/*!
 * \brief Calculate the harmonic mean of a fixed-size matrix.
 *
 * This is done by calculating the harmonic mean for each entry
 * individually. The result is stored the first argument.
 */
template <class Scalar, int m, int n>
void harmonicMeanMatrix(Dune::FieldMatrix<Scalar, m, n> &K,
                        const Dune::FieldMatrix<Scalar, m, n> &Ki,
                        const Dune::FieldMatrix<Scalar, m, n> &Kj)
{
    for (int rowIdx=0; rowIdx < m; rowIdx++){
        for (int colIdx=0; colIdx< n; colIdx++){
            if (Ki[rowIdx][colIdx] != Kj[rowIdx][colIdx]) {
                K[rowIdx][colIdx] = 
                    harmonicMean(Ki[rowIdx][colIdx],
                                 Kj[rowIdx][colIdx]);
            }
            else 
                K[rowIdx][colIdx] = Ki[rowIdx][colIdx];
        }
    }
}

}


#endif
