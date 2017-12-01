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
 * \brief This file contains free functions to evaluate the transmissibilities
 *        associated with flux evaluations across sub-control volume faces
 *        in the context of the cell-centered TPFA scheme.
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_COMPUTE_TRANSMISSIBILITY_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_COMPUTE_TRANSMISSIBILITY_HH

namespace Dumux
{

/*!
 * \ingroup Tpfa
 *
 * \brief Free function to evaluate the Tpfa transmissibility associated
 *        with the flux (in the form of flux = T*gradU) across a
 *        sub-control volume face stemming from a given sub-control
 *        volume with corresponding tensor T.
 *
 * \param scvf The sub-control volume face
 * \param scv The neighboring sub-control volume
 * \param K The tensor living in the neighboring scv
 * \param extrusionFactor The extrusion factor of the scv
 */
template< class SubControlVolumeFace, class SubControlVolume, class FieldScalar, int dimWorld >
FieldScalar computeTpfaTransmissibility(const SubControlVolumeFace& scvf,
                                        const SubControlVolume& scv,
                                        const Dune::FieldMatrix<FieldScalar, dimWorld, dimWorld>& T,
                                        typename SubControlVolume::Traits::Scalar extrusionFactor)
{
    using GlobalPosition = typename SubControlVolumeFace::Traits::GlobalPosition;
    GlobalPosition Knormal;
    T.mv(scvf.unitOuterNormal(), Knormal);

    auto distanceVector = scvf.ipGlobal();
    distanceVector -= scv.center();
    distanceVector /= distanceVector.two_norm2();

    return (Knormal*distanceVector) * extrusionFactor;
}

/*!
 * \ingroup Tpfa
 *
 * \brief Specialization of the above function for scalar T.
 *
 * \param scvf The sub-control volume face
 * \param scv The neighboring sub-control volume
 * \param t The scalar quantity living in the neighboring scv
 * \param extrusionFactor The extrusion factor of the scv
 */
template< class SubControlVolumeFace, class SubControlVolume, class FieldScalar >
FieldScalar computeTpfaTransmissibility(const SubControlVolumeFace& scvf,
                                        const SubControlVolume &scv,
                                        FieldScalar t,
                                        typename SubControlVolumeFace::Traits::Scalar extrusionFactor)
{
    auto distanceVector = scvf.ipGlobal();
    distanceVector -= scv.center();
    distanceVector /= distanceVector.two_norm2();

    return t * extrusionFactor * (distanceVector * scvf.unitOuterNormal());
}

} // end namespace Dumux

#endif
