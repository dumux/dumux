// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Hyperelastic
 * \brief Default implementation of the spatial params
 */
#ifndef DUMUX_SW_SOLIDMECHANICS_DEFAULT_HYPERELASTIC_SPATIAL_PARAMS_HH
#define DUMUX_SW_SOLIDMECHANICS_DEFAULT_HYPERELASTIC_SPATIAL_PARAMS_HH

#include <dumux/common/fvspatialparams.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class SandwichHyperelasticSpatialParams
: public FVSpatialParams<GridGeometry, Scalar, SandwichHyperelasticSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = FVSpatialParams<GridGeometry, Scalar, SandwichHyperelasticSpatialParams<GridGeometry, Scalar>>;
    using GlobalPosition = typename GridGeometry::SubControlVolume::GlobalPosition;
public:
    SandwichHyperelasticSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , E_(getParam<Scalar>("SpatialParams.YoungsModulus"))
    , nu_(getParam<Scalar>("SpatialParams.PoissonRatio"))
    , E2_(getParam<Scalar>("SpatialParams.YoungsModulus2"))
    , nu2_(getParam<Scalar>("SpatialParams.PoissonRatio2"))
    {
        mu_ = E_/(2*(1 + nu_));
        K_ = E_/(3*(1 - 2*nu_));
        lambda_ = nu_*E_/((1 + nu_)*(1-2*nu_));
        mu2_ = E2_/(2*(1 + nu2_));
        K2_ = E2_/(3*(1 - 2*nu2_));
        lambda2_ = nu2_*E2_/((1 + nu2_)*(1-2*nu2_));

        length_ = this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0];
    }

    Scalar shearModulus(const GlobalPosition& globalPos) const
    { return inMiddleLayer(globalPos) ? mu2_ : mu_; }

    Scalar bulkModulus(const GlobalPosition& globalPos) const
    { return inMiddleLayer(globalPos) ? K2_ : K_; }

    Scalar youngsModulus(const GlobalPosition& globalPos) const
    { return inMiddleLayer(globalPos) ? E2_ : E_; }

    Scalar poissonRatio(const GlobalPosition& globalPos) const
    { return inMiddleLayer(globalPos) ? nu2_ : nu_; }

    Scalar firstLameParameter(const GlobalPosition& globalPos) const
    { return inMiddleLayer(globalPos) ? lambda2_ : lambda_; }

private:
    Scalar E_, nu_, mu_, K_, lambda_;
    Scalar E2_, nu2_, mu2_, K2_, lambda2_;
    Scalar length_;

    bool inMiddleLayer(const GlobalPosition& globalPos) const
    {
        // Example condition: check if z-coordinate is within middle third of the domain height
        Scalar x = globalPos[0];
        return (x >= length_/3.0) && (x <= 2*length_/3.0);
    }
};

} // end namespace Dumux

#endif
