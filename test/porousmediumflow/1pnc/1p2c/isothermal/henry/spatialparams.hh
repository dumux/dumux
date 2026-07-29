// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Spatial parameters matching Fahs et al. (2016) Test Case 1 (Table 2).
 */
#ifndef DUMUX_TEST_HENRY_FAHS_SPATIAL_PARAMS_HH
#define DUMUX_TEST_HENRY_FAHS_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief Homogeneous, isotropic aquifer with the exact permeability and porosity
 *        given in Fahs et al. (2016), Table 2 (k = 1.0204e-9 m^2 is used directly,
 *        rather than re-derived from a hydraulic conductivity, to match their value
 *        exactly for this validation case).
 */
template<class GridGeometry, class Scalar>
class HenryFahsSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                             HenryFahsSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                           HenryFahsSpatialParams<GridGeometry, Scalar>>;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

public:
    using PermeabilityType = Scalar;

    HenryFahsSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

private:
    static constexpr Scalar permeability_ = 1.0204e-9; // [m^2], Fahs et al. 2016 Table 2
    static constexpr Scalar porosity_ = 0.35;
};

} // end namespace Dumux

#endif
