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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredFreeFlowDefaultFVGridGeometryTraits
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FREEFLOW_FV_GRID_GEOMETRY_TRAITS
#define DUMUX_DISCRETIZATION_STAGGERED_FREEFLOW_FV_GRID_GEOMETRY_TRAITS

#include <dumux/common/defaultmappertraits.hh>
#include <dumux/common/intersectionmapper.hh>
#include <dumux/discretization/cellcentered/subcontrolvolume.hh>
#include <dumux/discretization/staggered/fvelementgeometry.hh>
#include <dumux/discretization/staggered/freeflow/subcontrolvolumeface.hh>

#include "subcontrolvolumeface.hh"
#include "connectivitymap.hh"
#include "staggeredgeometryhelper.hh"

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Default traits for the finite volume grid geometry.
 */
template<class GridView, class MapperTraits = DefaultMapperTraits<GridView>>
struct StaggeredFreeFlowDefaultFVGridGeometryTraits
: public MapperTraits
{
    using SubControlVolume = CCSubControlVolume<GridView>;
    using SubControlVolumeFace = FreeFlowStaggeredSubControlVolumeFace<GridView>;
    using IntersectionMapper = ConformingGridIntersectionMapper<GridView>;
    using GeometryHelper = FreeFlowStaggeredGeometryHelper<GridView>;

    struct DofTypeIndices
    {
        using CellCenterIdx = Dune::index_constant<0>;
        using FaceIdx = Dune::index_constant<1>;
    };

    template<class FVGridGeometry>
    using ConnectivityMap = StaggeredFreeFlowConnectivityMap<FVGridGeometry>;

    template<class FVGridGeometry, bool cachingEnabled>
    using LocalView = StaggeredFVElementGeometry<FVGridGeometry, cachingEnabled>;
};

}

#endif
