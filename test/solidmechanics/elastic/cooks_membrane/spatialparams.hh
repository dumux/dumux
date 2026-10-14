// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsTests
 * \brief Spatial parameters of Cook's membrane: a homogeneous isotropic material given by Young's
 *        modulus and Poisson's ratio
 */
#ifndef DUMUX_TEST_ELASTIC_COOKS_MEMBRANE_SPATIAL_PARAMS_HH
#define DUMUX_TEST_ELASTIC_COOKS_MEMBRANE_SPATIAL_PARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/solidmechanics/elastic/lameparams.hh>
#include <dumux/solidmechanics/elastic/spatialparams.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class CooksMembraneSpatialParams
: public Experimental::ElasticSpatialParams<GridGeometry, Scalar, CooksMembraneSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = CooksMembraneSpatialParams<GridGeometry, Scalar>;
    using ParentType = Experimental::ElasticSpatialParams<GridGeometry, Scalar, ThisType>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

public:
    using LameParams = Dumux::LameParams<Scalar>;

    CooksMembraneSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        const auto E = getParam<Scalar>("SpatialParams.YoungsModulus");
        const auto nu = getParam<Scalar>("SpatialParams.PoissonRatio");
        lameParams_.setLambda(nu*E/((1.0 + nu)*(1.0 - 2.0*nu)));
        lameParams_.setMu(E/(2.0*(1.0 + nu)));
    }

    const LameParams& lameParamsAtPos(const GlobalPosition&) const
    { return lameParams_; }

private:
    LameParams lameParams_;
};

} // end namespace Dumux

#endif
