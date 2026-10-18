// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The spatial parameters of the vessel in the vessel tip benchmark.
 */
#ifndef DUMUX_TEST_EMBEDDED_1D3D_VESSEL_TIP_BLOODFLOW_SPATIALPARAMS_HH
#define DUMUX_TEST_EMBEDDED_1D3D_VESSEL_TIP_BLOODFLOW_SPATIALPARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief The spatial parameters of the vessel in the vessel tip benchmark.
 *
 * The vessel is a tube with cross-section \f$ \pi R^2 \f$ and conductivity
 * \f$ \kappa = \pi R^2 K \f$ (unit fluid viscosity).
 */
template<class GridGeometry, class Scalar>
class BloodFlowSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, BloodFlowSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = BloodFlowSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

public:
    using PermeabilityType = Scalar;

    BloodFlowSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        radius_ = getParam<Scalar>("SpatialParams.Radius");
        conductivity_ = getParam<Scalar>("SpatialParams.VesselConductivity");
    }

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return conductivity_/(M_PI*radius_*radius_); }

    Scalar radius(unsigned int eIdxGlobal) const
    { return radius_; }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    //! Extrudes the center-line to a tube with cross-section \f$ \pi R^2 \f$
    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    { return M_PI*radius_*radius_; }

private:
    Scalar radius_;
    Scalar conductivity_;
};

} // end namespace Dumux

#endif
