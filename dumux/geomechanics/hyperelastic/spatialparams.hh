// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Hyperelastic
 * \brief Default implementation of the spatial params
 */
#ifndef DUMUX_GEOMECHANICS_DEFAULT_HYPERELASTIC_SPATIAL_PARAMS_HH
#define DUMUX_GEOMECHANICS_DEFAULT_HYPERELASTIC_SPATIAL_PARAMS_HH

#include <dumux/common/fvspatialparams.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class DefaultHyperelasticSpatialParams
: public FVSpatialParams<GridGeometry, Scalar, DefaultHyperelasticSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = FVSpatialParams<GridGeometry, Scalar, DefaultHyperelasticSpatialParams<GridGeometry, Scalar>>;
public:
    DefaultHyperelasticSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , E_(getParam<Scalar>("Problem.YoungsModulus"))
    , nu_(getParam<Scalar>("Problem.PoissonRatio"))
    {
        mu_ = E_/(2*(1 + nu_));
        K_ = E_/(3*(1 - 2*nu_));
        lambda_ = nu_*E_/((1 + nu_)*(1-2*nu_));
    }

    Scalar shearModulus() const
    { return mu_; }

    Scalar bulkModulus() const
    { return K_; }

    Scalar youngsModulus() const
    { return E_; }

    Scalar poissonRatio() const
    { return nu_; }

    Scalar firstLameParameter() const
    { return lambda_; }

private:
    Scalar E_, nu_, mu_, K_, lambda_;
};

} // end namespace Dumux

#endif
