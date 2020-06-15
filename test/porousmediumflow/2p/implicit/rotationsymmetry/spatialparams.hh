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
 * \brief The spatial params 2p rotational symmetry test
 */
#ifndef DUMUX_TEST_TWOP_ROTATIONALSYMMETRY_SPATIAL_PARAMS_HH
#define DUMUX_TEST_TWOP_ROTATIONALSYMMETRY_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class TwoPRotationalSymmetrySpatialParams
: public FVSpatialParams<GridGeometry, Scalar, TwoPRotationalSymmetrySpatialParams<GridGeometry, Scalar>>
{
    using ThisType = TwoPRotationalSymmetrySpatialParams<GridGeometry, Scalar>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
public:
    using MaterialLaw = EffToAbsLaw<RegularizedBrooksCorey<Scalar>>;
    using MaterialLawParams = typename MaterialLaw::Params;
    using PermeabilityType = Scalar;

    TwoPRotationalSymmetrySpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // residual saturations
        materialLawParams_.setSwr(0.0);
        materialLawParams_.setSnr(0.0);

        // parameters for the material law
        materialLawParams_.setPe(0.0);
        materialLawParams_.setLambda(2.0);
    }

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return 1e-11; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    { return materialLawParams_; }

    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }

private:
    MaterialLawParams materialLawParams_;
};

} // end namespace Dumux

#endif
