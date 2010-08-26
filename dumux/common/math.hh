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
#ifndef DUMUX_MATH_HH
#define DUMUX_MATH_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <cmath>

namespace Dumux
{
/*!
 * \brief Calculate the harmonic mean of two scalar values.
 */
template <class Scalar>
Scalar harmonicMean(Scalar x, Scalar y)
{
    if (x*y <= 0)
        return 0;
    return (2*x*y)/(x + y);
}

/*!
 * \brief Calculate the geometric mean of two scalar values.
 */
template <class Scalar>
Scalar geometricMean(Scalar x, Scalar y)
{
    if (x*y <= 0)
        return 0;
    return std::sqrt(x*y)*((x < 0)?-1:1);
}

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

/*!
 * \brief Comparison of two position vectors
 *
 * Compares an current position vector with a reference vector, and returns true
 * if the position vector is larger.
 * "Larger" in this case means that all the entries of each spacial dimension are
 * larger compared to the reference vector.
 *
 * \param &pos Vektor holding the current Position that is to be checked
 * \return &smallerVec Reference vector, holding the minimum values for comparison.
 */
template <class Scalar, int dim>
bool isLarger(const Dune::FieldVector<Scalar, dim> &pos,
              const Dune::FieldVector<Scalar, dim> &smallerVec)
{
    for (int i=0; i < dim; i++)
    {
        if (pos[i]<= smallerVec[i])
        {
            return false;
        }
    }
    return true;
}

/*!
 * \brief Comparison of two position vectors
 *
 * Compares an current position vector with a reference vector, and returns true
 * if the position vector is smaller.
 * "Smaller" in this case means that all the entries of each spacial dimension are
 * smaller in comparison with the reference vector.
 *
 * \param &pos Vektor holding the current Position that is to be checked
 * \return &largerVec Reference vector, holding the maximum values for comparison.
 */
template <class Scalar, int dim>
bool isSmaller(const Dune::FieldVector<Scalar, dim> &pos,
               const Dune::FieldVector<Scalar, dim> &largerVec)
{
    for (int i=0; i < dim; i++)
    {
        if (pos[i]>= largerVec[i])
        {
            return false;
        }
    }
    return true;
}

/*!
 * \brief Comparison of three position vectors
 *
 * Compares an current position vector with two reference vector, and returns true
 * if the position vector lies in between them.
 * "Between" in this case means that all the entries of each spacial dimension are
 * smaller in comparison with the larger reference vector as well as larger campared
 * to the smaller reference.
 * This is comfortable to cheack weather the current position is located inside or
 * outside of a lense with different properties.
 *
 * \param &pos Vektor holding the current Position that is to be checked
 * \return &smallerVec Reference vector, holding the minimum values for comparison.
 * \return &largerVec Reference vector, holding the maximum values for comparison.
 */
template <class Scalar, int dim>
bool isBetween(const Dune::FieldVector<Scalar, dim> &pos,
              const Dune::FieldVector<Scalar, dim> &smallerVec,
              const Dune::FieldVector<Scalar, dim> &largerVec)
{
   if (isLarger(pos, smallerVec) && isSmaller(pos, largerVec))
       {
           return true;
       }
   else
       return false;
}


/*!
 * \brief Evaluates the Antoine equation used to calculate the vapour
 *        pressure of various liquids.
 *
 * See http://en.wikipedia.org/wiki/Antoine_equation
 */
template <class Scalar>
Scalar antoine(Scalar temperature,
               Scalar A,
               Scalar B,
               Scalar C)
{
    const Scalar ln10 = 2.3025850929940459;
    return std::exp(ln10*(A - B/(C + temperature)));
}

}


#endif
