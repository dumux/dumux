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
 * \brief Traits class to be used in conjunction with the CCMpfaFVGridGeometry.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FV_GRID_GEOMETRY_TRAITS_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FV_GRID_GEOMETRY_TRAITS_HH

#include <dumux/common/defaultmappertraits.hh>

#include <dumux/discretization/cellcentered/subcontrolvolume.hh>

#include <dumux/discretization/cellcentered/mpfa/connectivitymap.hh>
#include <dumux/discretization/cellcentered/mpfa/fvelementgeometry.hh>
#include <dumux/discretization/cellcentered/mpfa/subcontrolvolumeface.hh>
#include <dumux/discretization/cellcentered/mpfa/gridinteractionvolumeindexsets.hh>
#include <dumux/discretization/cellcentered/mpfa/helper.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Traits class to be used for the CCMpfaFVGridGeometry.
 *
 * \tparam GV the grid view type
 * \tparam NI the type used for node-local indexing
 * \tparam PIV the primary interaction volume type
 * \tparam SIV the secondary interaction volume type
 */
template<class GV, class NI, class PIV, class SIV>
struct CCMpfaFVGridGeometryTraits : public DefaultMapperTraits<GV>
{
    using SubControlVolume = CCSubControlVolume<GV>;
    using SubControlVolumeFace = CCMpfaSubControlVolumeFace<GV>;
    using NodalIndexSet = NI;
    //! State the maximum admissible element stencil size
    //! Per default, we use high values that are hopefully enough for all cases
    //! We assume simplex grids where stencils can get quite large but the size is unknown
    static constexpr int maxElementStencilSize = int(GV::dimension) == 3 ? 150 :
                                                 (int(GV::dimension)<int(GV::dimensionworld) ? 45 : 20);
    //! type definitions depending on the GridGeometry itself
    template< class FVGridGeom >
    using MpfaHelper = CCMpfaHelper< FVGridGeom >;
    template< class FVGridGeom >
    using ConnectivityMap = CCMpfaConnectivityMap<FVGridGeom, PIV::MpfaMethod>;
    template< class FVGridGeom >
    using GridIvIndexSets = CCMpfaGridInteractionVolumeIndexSets< FVGridGeom, NodalIndexSet, PIV, SIV >;
    template< class FVGridGeom, bool enableCache >
    using LocalView = CCMpfaFVElementGeometry<FVGridGeom, enableCache>;
};

} // end namespace Dumux

#endif
