// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
 * \param fvGeometry The element-centered control volume geometry
 * \param scvf The sub-control volume face
 * \param scv The neighboring sub-control volume
 * \param T The tensor living in the neighboring scv
 * \param extrusionFactor The extrusion factor of the scv
 */
template< class FVElementGeometry, class Tensor >
typename Tensor::field_type computeTpfaTransmissibility(const FVElementGeometry& fvGeometry,
                                                        const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                                        const typename FVElementGeometry::SubControlVolume& scv,
                                                        const Tensor& T,
                                                        typename FVElementGeometry::SubControlVolume::Traits::Scalar extrusionFactor)
{
    using GlobalPosition = typename FVElementGeometry::SubControlVolumeFace::Traits::GlobalPosition;
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
 * \param fvGeometry The element-centered control volume geometry
 * \param scvf The sub-control volume face
 * \param scv The neighboring sub-control volume
 * \param t The scalar quantity living in the neighboring scv
 * \param extrusionFactor The extrusion factor of the scv
 */
template< class FVElementGeometry,
          class Tensor,
          typename std::enable_if_t<Dune::IsNumber<Tensor>::value, int> = 0 >
Tensor computeTpfaTransmissibility(const FVElementGeometry& fvGeometry,
                                   const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                   const typename FVElementGeometry::SubControlVolume& scv,
                                   Tensor t,
                                   typename FVElementGeometry::SubControlVolumeFace::Traits::Scalar extrusionFactor)
{
    auto distanceVector = scvf.ipGlobal();
    distanceVector -= scv.center();
    distanceVector /= distanceVector.two_norm2();

    return t * extrusionFactor * (distanceVector * scvf.unitOuterNormal());
}

} // end namespace Dumux

#endif
