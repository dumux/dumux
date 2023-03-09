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

} // namespace Dumux

#endif
