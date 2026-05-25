// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_COMMON_EIGENVALUES_HH
#define DUMUX_COMMON_EIGENVALUES_HH

#include <cmath>
#include <limits>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

/*
 * \file
 * \ingroup Common
 * \brief Compute eigensystem of real symmetric 2x2 and 3x3 matrices
 * The algorithms in this header are adapted to Dumux from
 *
 * Numerical diagonalization of 3x3 matrcies
 * Copyright (C) 2006  Joachim Kopp, GNU LGPL v2.1 or later.
 *
 * If you are using this in a publication, consider citing:
 * Joachim Kopp, Numerical diagonalization of hermitian 3x3 matrices
 * arXiv.org preprint: physics/0610206
 * Int. J. Mod. Phys. C19 (2008) 523-548
 */

namespace Dumux::Detail {

/*!
 * \ingroup Common
 * \brief Compute eigensystem of real symmetric 2x2 matrix
 * Calculates the eigensystem of a real symmetric 2x2 matrix
 * [ A  B ] =  [ cosine  -sine ] [ eigenValue1   0  ] [  cosine  sine ]
 * [ B  C ]    [ sine   cosine ] [  0   eigenValue2 ] [ -sine  cosine ]
 * where eigenValue1 >= eigenValue2
 */
template<class Scalar>
inline void dsyev2(const Scalar A, const Scalar B, const Scalar C,
                   Scalar *eigenValue1, Scalar *eigenValue2,
                   Scalar *cosine, Scalar *sine)
{
    using std::sqrt; using std::abs;

    Scalar sm = A + C;
    Scalar df = A - C;
    Scalar rt = sqrt(df*df + 4.0*B*B);
    Scalar t;

    if (sm > 0.0)
    {
        *eigenValue1 = 0.5 * (sm + rt);
        t = 1.0/(*eigenValue1);
        *eigenValue2 = (A*t)*C - (B*t)*B;
    }
    else if (sm < 0.0)
    {
        *eigenValue2 = 0.5 * (sm - rt);
        t = 1.0/(*eigenValue2);
        *eigenValue1 = (A*t)*C - (B*t)*B;
    }
    // This case needs to be treated separately to avoid div by 0
    else
    {
        *eigenValue1 = 0.5 * rt;
        *eigenValue2 = -0.5 * rt;
    }

    // Calculate eigenvectors
    if (df > 0.0)
        *cosine = df + rt;
    else
        *cosine = df - rt;

    if (abs(*cosine) > 2.0*abs(B))
    {
        t = -2.0 * B / *cosine;
        *sine = 1.0 / sqrt(1.0 + t*t);
        *cosine = t * (*sine);
    }
    else if (abs(B) == 0.0)
    {
        *cosine = 1.0;
        *sine = 0.0;
    }
    else
    {
        t = -0.5 * (*cosine) / B;
        *cosine = 1.0 / sqrt(1.0 + t*t);
        *sine = t * (*cosine);
    }

    if (df > 0.0)
    {
        t = *cosine;
        *cosine = -(*sine);
        *sine = t;
    }
}

} // end namespace Dumux::Detail

namespace Dumux {

/*!
 * \ingroup Common
 * \brief Compute eigensystem of a real symmetric 2x2 matrix
 * \note eigenValue[0] >= eigenValue[1]
 */
template<class Scalar>
std::pair<Dune::FieldVector<Scalar, 2>, Dune::FieldMatrix<Scalar, 2, 2>>
eigenValuesVectors(const Dune::FieldMatrix<Scalar, 2, 2>& K)
{
    Dune::FieldVector<Scalar, 2> eigenValues;
    Dune::FieldMatrix<Scalar, 2, 2> eigenVectors;
    Dumux::Detail::dsyev2(K[0][0], K[1][0], K[1][1], &(eigenValues[0]), &(eigenValues[1]), &(eigenVectors[0][0]), &(eigenVectors[1][0]));
    eigenVectors[0][1] = -eigenVectors[1][0];
    eigenVectors[1][1] = eigenVectors[0][0];
    return std::make_pair(eigenValues, eigenVectors);
}

} // end namespace Dumux

#endif
