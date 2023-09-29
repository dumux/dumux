// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCMpfaDiscretization
 * \brief This file contains free functions to evaluate the transmissibilities
 *        associated with flux evaluations across sub-control volume faces
 *        in the context of cell-centered Mpfa schemes.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_COMPUTE_TRANSMISSIBILITY_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_COMPUTE_TRANSMISSIBILITY_HH

#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {

template<class T, class Enable = void>
struct TensorFieldType { using type = typename T::field_type; };

template<class T>
struct TensorFieldType<T, std::enable_if_t<Dune::IsNumber<T>::value>> { using type = T; };

} // namespace Detail
#endif // DOXYGEN

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Free function to evaluate the Mpfa transmissibility associated
 *        with the flux (in the form of flux = t*gradU) across a
 *        sub-control volume face stemming from a given sub-control
 *        volume with corresponding tensor t.
 *
 * \param fvGeometry The element-centered control volume geometry
 * \param scv The iv-local sub-control volume
 * \param scvf The grid sub-control volume face
 * \param t The tensor/Scalar living in the scv
 * \param extrusionFactor The extrusion factor of the scv
 */
template<class EG, class IVSubControlVolume, class Tensor>
auto computeMpfaTransmissibility(const EG& fvGeometry,
                                 const IVSubControlVolume& scv,
                                 const typename EG::SubControlVolumeFace& scvf,
                                 const Tensor& t,
                                 typename IVSubControlVolume::ctype extrusionFactor)
{
    Dune::FieldVector<typename Detail::TensorFieldType<Tensor>::type, IVSubControlVolume::myDimension> wijk;
    for (unsigned int dir = 0; dir < IVSubControlVolume::myDimension; ++dir)
        wijk[dir] = vtmv(scvf.unitOuterNormal(), t, scv.nu(dir));

    using Extrusion = Extrusion_t<typename EG::GridGeometry>;
    wijk *= Extrusion::area(fvGeometry, scvf)*extrusionFactor;
    wijk /= scv.detX();

    return wijk;
}

} // end namespace Dumux

#endif
