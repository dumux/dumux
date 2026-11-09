// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Constant spatial parameters with dispersivities for the Henry problem \cite Henry1964
 */
#ifndef DUMUX_TEST_HENRY_FAHS_SPATIALPARAMS_HH
#define DUMUX_TEST_HENRY_FAHS_SPATIALPARAMS_HH

#include <array>
#include <memory>

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/fvspatialparams1pconstant.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief Constant spatial parameters with the longitudinal and transverse dispersivities
 *        of Scheidegger's dispersion tensor
 */
template<class GridGeometry, class Scalar>
class HenryFahsSpatialParams
: public FVPorousMediumFlowSpatialParamsOnePConstant<GridGeometry, Scalar>
{
    using ParentType = FVPorousMediumFlowSpatialParamsOnePConstant<GridGeometry, Scalar>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Geometry::GlobalCoordinate;

public:
    HenryFahsSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , dispersivities_{getParam<Scalar>("SpatialParams.LongitudinalDispersivity"),
                      getParam<Scalar>("SpatialParams.TransverseDispersivity")}
    {}

    //! The longitudinal and transverse dispersivity \f$[m]\f$
    std::array<Scalar, 2> dispersionAlphas(const GlobalPosition& globalPos, int phaseIdx = 0, int compIdx = 0) const
    { return dispersivities_; }

private:
    std::array<Scalar, 2> dispersivities_;
};

} // end namespace Dumux

#endif
