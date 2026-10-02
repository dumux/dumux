// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTests
 * \brief Spatial parameters for the pressure step of the IMPES Buckley-Leverett test.
 *
 * The pressure equation of the IMPES scheme
 * \f[ -\nabla\cdot\left(K \lambda_t(S_w^n) \nabla p\right) = 0 \f]
 * is solved with the single-phase model. To this end, the single-phase permeability
 * is replaced by an effective permeability \f$ K \lambda_t \mu_\text{ref} \f$, where
 * \f$ \mu_\text{ref} \f$ is the viscosity of the fluid used by the single-phase model.
 * The factor \f$ \lambda_t \mu_\text{ref} \f$ is computed from the old saturation and
 * set from outside via setRelativeTotalMobility().
 *
 * The intrinsic permeability, porosity and temperature are taken from the two-phase
 * spatial parameters, so both models use the same medium. This class is kept separate
 * on purpose: its permeability() is an effective permeability including the total
 * mobility, which must not be used by the two-phase model (which multiplies the
 * intrinsic permeability with the phase mobilities itself).
 *
 * \tparam TwoPSpatialParams the spatial parameters of the two-phase problem
 */
#ifndef DUMUX_TEST_TWOP_BUCKLEYLEVERETT_SPATIALPARAMS_IMPES_HH
#define DUMUX_TEST_TWOP_BUCKLEYLEVERETT_SPATIALPARAMS_IMPES_HH

#include <memory>
#include <vector>

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

template<class GridGeometry, class Scalar, class TwoPSpatialParams>
class BuckleyLeverettImpesPressureSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, BuckleyLeverettImpesPressureSpatialParams<GridGeometry, Scalar, TwoPSpatialParams>>
{
    using ThisType = BuckleyLeverettImpesPressureSpatialParams<GridGeometry, Scalar, TwoPSpatialParams>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = typename TwoPSpatialParams::PermeabilityType;

    BuckleyLeverettImpesPressureSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , twoPSpatialParams_(gridGeometry)
    , relativeTotalMobility_(gridGeometry->numDofs(), 1.0)
    {}

    //! The effective permeability \f$ K \lambda_t \mu_\text{ref} \f$
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        PermeabilityType k = twoPSpatialParams_.permeability(element, scv, elemSol);
        k *= relativeTotalMobility_[scv.dofIndex()];
        return k;
    }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return twoPSpatialParams_.porosityAtPos(globalPos); }

    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return twoPSpatialParams_.temperatureAtPos(globalPos); }

    //! Set the factor \f$ \lambda_t \mu_\text{ref} \f$ for each cell
    void setRelativeTotalMobility(const std::vector<Scalar>& relativeTotalMobility)
    { relativeTotalMobility_ = relativeTotalMobility; }

    const std::vector<Scalar>& relativeTotalMobility() const
    { return relativeTotalMobility_; }

private:
    TwoPSpatialParams twoPSpatialParams_;
    std::vector<Scalar> relativeTotalMobility_;
};

} // end namespace Dumux

#endif
