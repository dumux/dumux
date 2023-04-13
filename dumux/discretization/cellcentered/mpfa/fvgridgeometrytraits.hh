// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
