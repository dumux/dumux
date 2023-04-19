// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \ingroup SpatialParameters
 * \brief Spatial parameters for an isothermal 1p pore-network model
 */
#ifndef DUMUX_COMPOSITIONAL_PNM_SPATIAL_PARAMS_1P_HH
#define DUMUX_COMPOSITIONAL_PNM_SPATIAL_PARAMS_1P_HH

#include <dumux/porenetwork/common/spatialparams.hh>

namespace Dumux::PoreNetwork {

template<class GridGeometry, class Scalar>
class CompositionalSpatialParams
: public SpatialParams<GridGeometry, Scalar, CompositionalSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = SpatialParams<GridGeometry, Scalar, CompositionalSpatialParams<GridGeometry, Scalar>>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
public:
    using PermeabilityType = Scalar;
    using ParentType::ParentType;

    CompositionalSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    { }

    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 273.15 + 10.0; }
};
} // end namespace Dumux::PoreNetwork

#endif
