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
 * \ingroup CCTpfaDiscretization
 * \brief Free functions to evaluate the transmissibilities
 *        associated with flux evaluations across sub-control volume faces
 *        in the context of the cell-centered TPFA scheme.
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_COMPUTE_TRANSMISSIBILITY_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_COMPUTE_TRANSMISSIBILITY_HH

#include <dune/common/typetraits.hh>
#include <dune/common/fmatrix.hh>

namespace Dumux {

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Free function to evaluate the Tpfa transmissibility
 *        associated with the flux (in the form of flux = T*gradU) across a
 *        sub-control volume face stemming from a given sub-control
 *        volume with corresponding tensor T.
 *
 * \param scvf The sub-control volume face
 * \param scv The neighboring sub-control volume
 * \param T The tensor living in the neighboring scv
 * \param extrusionFactor The extrusion factor of the scv
 */
template< class SubControlVolumeFace, class SubControlVolume, class Tensor >
typename Tensor::field_type computeTpfaTransmissibility(const SubControlVolumeFace& scvf,
                                                        const SubControlVolume& scv,
                                                        const Tensor& T,
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
 * \ingroup CCTpfaDiscretization
 * \brief Free function to evaluate the Tpfa transmissibility
 *        associated with the flux (in the form of flux = T*gradU) across a
 *        sub-control volume face stemming from a given sub-control
 *        volume for the case where T is just a scalar
 *
 * \param scvf The sub-control volume face
 * \param scv The neighboring sub-control volume
 * \param t The scalar quantity living in the neighboring scv
 * \param extrusionFactor The extrusion factor of the scv
 */
template< class SubControlVolumeFace,
          class SubControlVolume,
          class Tensor,
          typename std::enable_if_t<Dune::IsNumber<Tensor>::value, int> = 0 >
Tensor computeTpfaTransmissibility(const SubControlVolumeFace& scvf,
                                   const SubControlVolume &scv,
                                   Tensor t,
                                   typename SubControlVolumeFace::Traits::Scalar extrusionFactor)
{
    auto distanceVector = scvf.ipGlobal();
    distanceVector -= scv.center();
    distanceVector /= distanceVector.two_norm2();

    return t * extrusionFactor * (distanceVector * scvf.unitOuterNormal());
}

} // end namespace Dumux

#endif
