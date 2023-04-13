// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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

#include "subcontrolvolumeface.hh"
#include "connectivitymap.hh"
#include "staggeredgeometryhelper.hh"

namespace Dumux {

template<class GridView>
class FreeflowStaggeredSCV
{
    using ctype = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<ctype, GridView::dimensionworld>;
public:
    struct Traits
    {
        using Scalar = ctype;
        using Geometry = typename GridView::template Codim<0>::Geometry;
    };

    FreeflowStaggeredSCV(const GlobalPosition& center, const ctype volume)
    : center_(center), volume_(volume) {}

    const GlobalPosition& center() const
    { return center_; }

    ctype volume() const
    { return volume_; }
private:
    GlobalPosition center_;
    ctype volume_;
};

template<class GridView>
class FreeflowStaggeredSCVF
{
    using ctype = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<ctype, GridView::dimensionworld>;
public:
    struct Traits
    {
        using Scalar = ctype;
        using Geometry = typename GridView::template Codim<1>::Geometry;
    };

    FreeflowStaggeredSCVF(const GlobalPosition& center, const ctype area)
    : center_(center), area_(area) {}

    const GlobalPosition& center() const
    { return center_; }

    ctype area() const
    { return area_; }
private:
    GlobalPosition center_;
    ctype area_;
};


/*!
 * \ingroup StaggeredDiscretization
 * \brief Default traits for the finite volume grid geometry.
 */
template<class GridView, int upwOrder, class MapperTraits = DefaultMapperTraits<GridView>>
struct StaggeredFreeFlowDefaultFVGridGeometryTraits
: public MapperTraits
{
    using SubControlVolume = CCSubControlVolume<GridView>;
    using SubControlVolumeFace = FreeFlowStaggeredSubControlVolumeFace<GridView, upwOrder>;
    using IntersectionMapper = ConformingGridIntersectionMapper<GridView>;
    using GeometryHelper = FreeFlowStaggeredGeometryHelper<GridView, upwOrder>;
    static constexpr int upwindSchemeOrder = upwOrder;

    struct DofTypeIndices
    {
        using FaceIdx = Dune::index_constant<0>;
        using CellCenterIdx = Dune::index_constant<1>;
    };

    template<class GridGeometry>
    using ConnectivityMap = StaggeredFreeFlowConnectivityMap<GridGeometry>;

    template<class GridGeometry, bool cachingEnabled>
    using LocalView = StaggeredFVElementGeometry<GridGeometry, cachingEnabled>;

    struct PublicTraits
    {
        using CellSubControlVolume = SubControlVolume;
        using CellSubControlVolumeFace = SubControlVolumeFace;
        using FaceSubControlVolume = FreeflowStaggeredSCV<GridView>;
        using FaceLateralSubControlVolumeFace = FreeflowStaggeredSCVF<GridView>;
        using FaceFrontalSubControlVolumeFace = FreeflowStaggeredSCVF<GridView>;
    };
};

} //end namespace Dumux

#endif
