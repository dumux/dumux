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
 *
 * \brief Definition of the spatial parameters for the 1p2c
 *       salt instability problem
 */
#ifndef DUMUX_1P2C_SPATIAL_PARAMETERS_HH
#define DUMUX_1P2C_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux {

/*!
 *
 * \brief Definition of the spatial parameters for the 1p2c
 *        salt instability problem
 */
template<class FVGridGeometry, class Scalar>
class SaltInstability1p2cSpatialParams :
public FVSpatialParams<FVGridGeometry, Scalar,
                       SaltInstability1p2cSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar,
                                       SaltInstability1p2cSpatialParams<FVGridGeometry, Scalar>>;

    using GlobalPosition = typename SubControlVolume::GlobalPosition;

public:
    using PermeabilityType = Scalar;

    SaltInstability1p2cSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        // intrinsic permeabilities
        perm_ = 1.019368e-8;
        porosity_ = 0.35;
    }

    /*!
     * \brief Define the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        return perm_;
    }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    double porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_;
    }


private:
    Scalar perm_;
    Scalar porosity_;
};

} // end namespace Dumux

#endif // DUMUX_HENRY1P2C_SPATIAL_PARAMETERS_HH
