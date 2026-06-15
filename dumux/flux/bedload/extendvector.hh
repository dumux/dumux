// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadFlux
 * \copydoc Dumux::Sediment::modifyVectorByAngle
 *
 */
#ifndef DUMUX_FLUX_MODIFY_VECTOR_HH
#define DUMUX_FLUX_MODIFY_VECTOR_HH

namespace Dumux {
namespace Sediment {
/*!
 * \ingroup BedloadFlux
 * \brief Return a extended vector
 *
 * The vector is extended by a perpendicular component. This perpendicular component is not stated explicitly.
 * Instead the tangent between the perpendicular component and the original vector is used to quantify the modidication.
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \param vector The vector which shall be modified. The size must be (2 x 1).
 * \param tan_angle Tangent between the perpendicular component and the original vector.
 *
 * \return An extended copy of the given vector
 */
template<class Scalar>
const Dune::FieldVector<Scalar, 2> getExtendedVector(Dune::FieldVector<Scalar, 2> vector, Scalar tan_angle)
{
    Dune::FieldVector<Scalar, 2> vector_modified(0.0);

    vector_modified[0] = vector[0] - tan_angle * vector[1];
    vector_modified[1] = vector[1] + tan_angle * vector[0];

    return vector_modified;
}
} // end namespace Sediment
} // end namespace Dumux

#endif
