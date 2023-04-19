// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup SolidEnergyTests
 * \brief Definition of the spatial parameters for the solid energy test
 */
#ifndef DUMUX_TEST_SOLIDENERGY_SPATIAL_PARAMS_HH
#define DUMUX_TEST_SOLIDENERGY_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief Definition of the spatial parameters for the solid energy test
 */
template<class GridGeometry, class Scalar>
class SolidEnergySpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                         SolidEnergySpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                                       SolidEnergySpatialParams<GridGeometry, Scalar>>;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    SolidEnergySpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry) {}

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }
};

} // end namespace Dumux

#endif
