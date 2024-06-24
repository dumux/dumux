// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
 * \note We do a harmonic average of the part normal to the interface (alpha nn^T) and
 *       an arithmetic average of the tangential part (T - alpha nn^T).
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
    const Scalar alpha1 = vtmv(normal, T1, normal);
    const Scalar alpha2 = vtmv(normal, T2, normal);

    const auto normalProjection = dyadicProduct(normal, normal);

    auto T = (Dumux::harmonicMean(alpha1, alpha2) - 0.5*(alpha1+alpha2))*normalProjection
           + 0.5*(T1 + T2);

    return T;
}

} // end namespace Dumux

#endif
