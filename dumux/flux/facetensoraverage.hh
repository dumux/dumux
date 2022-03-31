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
 * \ingroup Flux
 * \brief A free function to average a Tensor at an interface.
 */
#ifndef DUMUX_FACE_TENSOR_AVERAGE_HH
#define DUMUX_FACE_TENSOR_AVERAGE_HH

#include <dune/common/fmatrix.hh>
#include <dumux/common/math.hh>

namespace Dumux {

/*!
 * \brief Average of a discontinuous scalar field at discontinuity interface
 *        (for compatibility reasons with the function below)
 * \return the harmonic average of the scalars
 * \param T1 first scalar parameter
 * \param T2 second scalar parameter
 * \param normal The unit normal vector of the interface
 */
template <class Scalar, int dim>
Scalar faceTensorAverage(const Scalar T1,
                         const Scalar T2,
                         const Dune::FieldVector<Scalar, dim>& normal)
{ return Dumux::harmonicMean(T1, T2); }

/*!
 * \brief Average of a discontinuous tensorial field at discontinuity interface
 * \note We do a harmonic average of the part normal to the interface (alpha*I) and
 *       an arithmetic average of the tangential part (T - alpha*I).
 * \return the averaged tensor
 * \param T1 first tensor
 * \param T2 second tensor
 * \param normal The unit normal vector of the interface
 */
template <class Scalar, int dim>
Dune::FieldMatrix<Scalar, dim> faceTensorAverage(const Dune::FieldMatrix<Scalar, dim>& T1,
                                                 const Dune::FieldMatrix<Scalar, dim>& T2,
                                                 const Dune::FieldVector<Scalar, dim>& normal)
{
    // determine nT*k*n
    Dune::FieldVector<Scalar, dim> tmp;
    Dune::FieldVector<Scalar, dim> tmp2;
    T1.mv(normal, tmp);
    T2.mv(normal, tmp2);
    const Scalar alpha1 = tmp*normal;
    const Scalar alpha2 = tmp2*normal;

    const Scalar alphaHarmonic = Dumux::harmonicMean(alpha1, alpha2);
    const Scalar alphaAverage = 0.5*(alpha1 + alpha2);

    Dune::FieldMatrix<Scalar, dim> T(0.0);
    for (int i = 0; i < dim; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            T[i][j] += 0.5*(T1[i][j] + T2[i][j]);
            if (i == j)
                T[i][j] += alphaHarmonic - alphaAverage;
        }
    }

    return T;
}

} // end namespace Dumux

#endif
