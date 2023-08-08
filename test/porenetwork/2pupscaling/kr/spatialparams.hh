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
 * \brief Spatial parameters for an isothermal 2p pore-network model
 */
#ifndef DUMUX_PNM_2P_DRAINAGE_SPATIAL_PARAMS_HH
#define DUMUX_PNM_2P_DRAINAGE_SPATIAL_PARAMS_HH

#include <dumux/porenetwork/2p/spatialparams.hh>

namespace Dumux::PoreNetwork {

template<class GridGeometry, class Scalar, class MaterialLawT>
class TwoPDrainageSpatialParams
: public TwoPSpatialParams<GridGeometry, Scalar, MaterialLawT, TwoPDrainageSpatialParams<GridGeometry, Scalar, MaterialLawT>>
{
    using ParentType = TwoPSpatialParams<GridGeometry, Scalar, MaterialLawT, TwoPDrainageSpatialParams<GridGeometry, Scalar, MaterialLawT>>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
public:
    using PermeabilityType = Scalar;
    using ParentType::ParentType;

    TwoPDrainageSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        temperature_ = getParam<Scalar>("SpatialParams.Temperature", 283.15);
        gamma_ = getParam<Scalar>("SpatialParams.SurfaceTension", 0.0725); // default to surface tension of water/air
        theta_ = getParam<Scalar>("SpatialParams.ContactAngle", 0.0);
    }

    //! \brief Function for defining the temperature
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return temperature_; }

    //! \brief Function for defining the wetting phase
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }


    //! \brief Function for defining the contact angle
    int contactAngleAtPos(const GlobalPosition& globalPos) const
    { return theta_; }

    //! \brief Function for defining the surface tension
    Scalar surfaceTensionAtPos(const GlobalPosition& globalPos) const
    { return gamma_; }

private:
    Scalar temperature_;
    Scalar gamma_;
    Scalar theta_;
};
} // end namespace Dumux::PoreNetwork

#endif
