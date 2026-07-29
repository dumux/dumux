// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Spatial parameters matching Fahs et al. (2016), Table 2.
 */
#ifndef DUMUX_TEST_HENRY_FAHS_SPATIAL_PARAMS_HH
#define DUMUX_TEST_HENRY_FAHS_SPATIAL_PARAMS_HH

#include <array>

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief Homogeneous, isotropic aquifer with the exact permeability and porosity
 *        given in Fahs et al. (2016), Table 2 (k = 1.0204e-9 m^2 is used directly,
 *        rather than re-derived from a hydraulic conductivity, to match their value
 *        exactly). Longitudinal/transverse dispersivities (alphaL, alphaT) are read
 *        from Problem.AlphaL/Problem.AlphaT: 0/0 for "Test Case 1" (purely diffusive,
 *        dispersion disabled entirely -- see properties.hh), nonzero for Test Cases
 *        2/3 where velocity-dependent (Scheidegger) dispersion is enabled.
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
    , alphaL_(getParam<Scalar>("Problem.AlphaL"))
    , alphaT_(getParam<Scalar>("Problem.AlphaT"))
    {}

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    //! Longitudinal/transverse dispersivities [m] for Scheidegger's dispersion tensor
    //! (only used when EnableCompositionalDispersion is set, see properties.hh)
    std::array<Scalar, 2> dispersionAlphas(const GlobalPosition& globalPos, int phaseIdx = 0, int compIdx = 0) const
    { return {alphaL_, alphaT_}; }

private:
    static constexpr Scalar permeability_ = 1.0204e-9; // [m^2], Fahs et al. 2016 Table 2
    static constexpr Scalar porosity_ = 0.35;
    Scalar alphaL_;
    Scalar alphaT_;
};

} // end namespace Dumux

#endif
