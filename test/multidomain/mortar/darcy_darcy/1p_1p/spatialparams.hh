// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \ingroup MultiDomainTests
 * \brief The spatial parameters for the single-phase darcy-darcy mortar-coupling test
 */
#ifndef DUMUX_MORTAR_DARCY_ONEP_TEST_SPATIAL_PARAMS_HH
#define DUMUX_MORTAR_DARCY_ONEP_TEST_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

#ifndef PERMEABILITY_JUMP
#define PERMEABILITY_JUMP 0
#endif


/*!
 * \ingroup OnePTests
 * \ingroup MultiDomainTests
 * \brief The spatial parameters for the single-phase darcy-darcy mortar-coupling test
 */
template<class GridGeometry, class Scalar>
class DarcySpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                             DarcySpatialParams<GridGeometry, Scalar>>
{
    using ThisType = DarcySpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;
    DarcySpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , permeability_{getParam<Scalar>("SpatialParams.Permeability")}
    {}

    Scalar permeabilityAtPos(const GlobalPosition& globalPos) const
    {
#if PERMEABILITY_JUMP
        // each subdomain lies wholly on one side of the interface, so this is constant per subdomain
        return globalPos[0] < 0.5 ? 1.0 : 2.0;
#else
        return permeability_;
#endif
    }
    Scalar porosityAtPos(const GlobalPosition& globalPos) const { return 1.0; }
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const { return 283.15; }
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const { return 1.0; }

private:
    Scalar permeability_;
};

} // end namespace Dumux

#endif
