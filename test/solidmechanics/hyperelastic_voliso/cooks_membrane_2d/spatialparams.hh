// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsTests
 * \brief Spatial parameters of Cook's membrane
 */
#ifndef DUMUX_HYPERELASTIC_COOKS_MEMBRANE_SPATIAL_PARAMS_HH
#define DUMUX_HYPERELASTIC_COOKS_MEMBRANE_SPATIAL_PARAMS_HH

#include <dumux/common/spatialparams.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup GeomechanicsTests
 * \brief Lamé parameters of Cook's membrane, read from the input
 *
 * Lamé's first parameter can be changed between solves, which the tests use to run a sequence of
 * increasingly incompressible materials.
 */
template<class GridGeometry, class Scalar>
class CooksMembraneHyperelasticSpatialParams
: public Experimental::SpatialParams<GridGeometry, Scalar,
                                     CooksMembraneHyperelasticSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = Experimental::SpatialParams<GridGeometry, Scalar,
                                                   CooksMembraneHyperelasticSpatialParams<GridGeometry, Scalar>>;
public:
    CooksMembraneHyperelasticSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , lambda_(getParam<Scalar>("SpatialParams.Lambda"))
    , mu_(getParam<Scalar>("SpatialParams.Mu"))
    {}

    Scalar shearModulus() const { return mu_; }
    Scalar firstLameParameter() const { return lambda_; }
    Scalar bulkModulus() const { return lambda_ + 2.0*mu_/3.0; }

    void setLambda(Scalar lambda) { lambda_ = lambda; }

private:
    Scalar lambda_, mu_;
};

} // end namespace Dumux

#endif
