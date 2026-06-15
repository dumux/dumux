// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadFlux
 * \copydoc Dumux::Sediment::rotateVector
 *
 */
#ifndef DUMUX_FLUX_ROTATE_VECTOR_HH
#define DUMUX_FLUX_ROTATE_VECTOR_HH

namespace Dumux {
namespace Sediment {
/*!
 * \ingroup BedloadFlux
 * \brief Return a rotated vector
 *
 * The vector is rotated counterclockwise by an angle.
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \param vector The vector which shall be rotated. The size must be (2 x 1).
 * \param angle Angle (radian measure) by which the vector will be rotated.
 *
 * \return A rotated copy of the given vector
 */
template<class Scalar>
const Dune::FieldVector<Scalar, 2> getRotatedVector(Dune::FieldVector<Scalar, 2> vector, Scalar angle)
{
    using std::sin;
    using std::cos;
    Dune::FieldVector<Scalar, 2> vector_rotated(0.0);

    vector_rotated[0] = cos(angle) * vector[0] - sin(angle) * vector[1];
    vector_rotated[1] = sin(angle) * vector[0] + cos(angle) * vector[1];

    return vector_rotated;
}
} // end namespace Sediment
} // end namespace Dumux

#endif
