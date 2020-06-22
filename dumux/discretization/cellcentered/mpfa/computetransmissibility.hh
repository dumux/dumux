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

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Free function to evaluate the Mpfa transmissibility associated
 *        with the flux (in the form of flux = t*gradU) across a
 *        sub-control volume face stemming from a given sub-control
 *        volume with corresponding tensor t.
 *
 * \param scv The iv-local sub-control volume
 * \param scvf The grid sub-control volume face
 * \param t The tensor living in the scv
 * \param extrusionFactor The extrusion factor of the scv
 */
template<class EG, class IVSubControlVolume, class Tensor>
Dune::FieldVector<typename Tensor::field_type, IVSubControlVolume::myDimension>
computeMpfaTransmissibility(const IVSubControlVolume& scv,
                            const typename EG::SubControlVolumeFace& scvf,
                            const Tensor& t,
                            typename IVSubControlVolume::ctype extrusionFactor)
{
    Dune::FieldVector<typename Tensor::field_type, IVSubControlVolume::myDimension> wijk;
    for (unsigned int dir = 0; dir < IVSubControlVolume::myDimension; ++dir)
        wijk[dir] = vtmv(scvf.unitOuterNormal(), t, scv.nu(dir));

    using Extrusion = Extrusion_t<typename EG::GridGeometry>;
    wijk *= Extrusion::area(scvf)*extrusionFactor;
    wijk /= scv.detX();

    return wijk;
}

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Free function to evaluate the Mpfa transmissibility associated
 *        with the flux (in the form of flux = t*gradU) across a
 *        sub-control volume face stemming from a given sub-control
 *        volume with corresponding tensor t, where t is a scalar.
 *
 * \param scv The iv-local sub-control volume
 * \param scvf The grid sub-control volume face
 * \param t the scalar quantity living in the scv
 * \param extrusionFactor The extrusion factor of the scv
 */
template< class EG, class IVSubControlVolume, class Tensor, std::enable_if_t< Dune::IsNumber<Tensor>::value, int > = 1 >
Dune::FieldVector<Tensor, IVSubControlVolume::myDimension>
computeMpfaTransmissibility(const IVSubControlVolume& scv,
                            const typename EG::SubControlVolumeFace& scvf,
                            const Tensor& t,
                            typename IVSubControlVolume::ctype extrusionFactor)
{
    Dune::FieldVector<Tensor, IVSubControlVolume::myDimension> wijk;
    for (unsigned int dir = 0; dir < IVSubControlVolume::myDimension; ++dir)
        wijk[dir] = vtmv(scvf.unitOuterNormal(), t, scv.nu(dir));

    using Extrusion = Extrusion_t<typename EG::GridGeometry>;
    wijk *= Extrusion::area(scvf)*extrusionFactor;
    wijk /= scv.detX();

    return wijk;
}

} // end namespace Dumux

#endif
