// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMSolidEnergyModel
 * \ingroup SpatialParameters
 * \brief The spatial parameters for solid-energy models in pore networks.
 */
#ifndef DUMUX_PNM_SOLID_ENERGY_SPATIAL_PARAMS_HH
#define DUMUX_PNM_SOLID_ENERGY_SPATIAL_PARAMS_HH

#include <dumux/porenetwork/common/spatialparams.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMSolidEnergyModel
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters for pore network solid models.
 */
template<class FVGridGeometry, class Scalar>
class SolidEnergySpatialParams
: public PoreNetwork::SpatialParams<FVGridGeometry, Scalar,
                                        SolidEnergySpatialParams<FVGridGeometry, Scalar>>
{
    using ParentType = PoreNetwork::SpatialParams<
        FVGridGeometry, Scalar, SolidEnergySpatialParams<FVGridGeometry, Scalar>
    >;

    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

public:
    using PermeabilityType = Scalar;
    using ParentType::ParentType;

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.0; }

};

} // end namespace Dumux::PoreNetwork

#endif
