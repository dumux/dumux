// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_WOOD_TEST_SPATIAL_PARAMS_HH
#define DUMUX_WOOD_TEST_SPATIAL_PARAMS_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/fvspatialparams.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \brief Spatial parameters for the orthotropic wood test (transverse \f$R\f$-\f$T\f$ cross-section).
 *
 * The 2D plane is perpendicular to the longitudinal grain (\f$L\f$ is out-of-plane).
 * \f$R\f$ = radial direction (from pith outward), \f$T\f$ = tangential (along annual rings).
 * The local frame \f$(R,T)\f$ varies with position: at position \f$\mathbf{p}\f$
 * \f[
 *     \mathbf{e}_R(\mathbf{p}) = \frac{\mathbf{p} - \mathbf{p}_\mathrm{pith}}{\|\mathbf{p} - \mathbf{p}_\mathrm{pith}\|},
 *     \qquad \mathbf{e}_T = \mathbf{R}_{90}\, \mathbf{e}_R.
 * \f]
 * Stiffness coefficients are stored in the local frame; rotation to the global
 * frame is applied in the problem's constitutive law.
 */
template<class GridGeometry, class Scalar>
class WoodSpatialParams
: public FVSpatialParams<GridGeometry, Scalar, WoodSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = WoodSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

public:
    using RotationMatrix = Dune::FieldMatrix<Scalar, 2, 2>;

    WoodSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , eR_(getParam<Scalar>("SpatialParams.YoungsModulusR"))
    , eT_(getParam<Scalar>("SpatialParams.YoungsModulusT"))
    , nuRT_(getParam<Scalar>("SpatialParams.PoissonRatioRT"))
    , gRT_(getParam<Scalar>("SpatialParams.ShearModulusRT"))
    , alphaR_(getParam<Scalar>("SpatialParams.ShrinkageR"))
    , alphaT_(getParam<Scalar>("SpatialParams.ShrinkageT"))
    , mRef_(getParam<Scalar>("SpatialParams.ReferenceMoisture"))
    , dR_(getParam<Scalar>("SpatialParams.MoistureDiffusivityR"))
    , dT_(getParam<Scalar>("SpatialParams.MoistureDiffusivityT"))
    , hM_(getParam<Scalar>("SpatialParams.MassTransferCoefficient"))
    , mAir_(getParam<Scalar>("SpatialParams.AirMoisture"))
    , mInit_(getParam<Scalar>("SpatialParams.InitialMoisture"))
    , pith_({getParam<Scalar>("SpatialParams.PithX"),
             getParam<Scalar>("SpatialParams.PithY")})
    {
        const Scalar nuTR = nuRT_ * eT_ / eR_;        // Maxwell-Betti
        const Scalar denom = 1.0 - nuRT_ * nuTR;
        c11_ = eR_ / denom;                           // local axis 1 = R
        c22_ = eT_ / denom;                           // local axis 2 = T
        c12_ = nuRT_ * eT_ / denom;                   // = nuTR * eR_ / denom
        c66_ = gRT_;
    }

    Scalar C11() const { return c11_; }
    Scalar C22() const { return c22_; }
    Scalar C12() const { return c12_; }
    Scalar C66() const { return c66_; }

    Scalar shrinkageR() const { return alphaR_; }
    Scalar shrinkageT() const { return alphaT_; }
    Scalar referenceMoisture() const { return mRef_; }

    // anisotropic moisture diffusivity in the global frame at globalPos
    // D_global(x) = Q(x) · diag(D_R, D_T) · Q(x)^T
    RotationMatrix moistureDiffusivity(const GlobalPosition& globalPos) const
    {
        const auto Q = rotation(globalPos);
        RotationMatrix D{{dR_*Q[0][0]*Q[0][0] + dT_*Q[0][1]*Q[0][1],
                          dR_*Q[0][0]*Q[1][0] + dT_*Q[0][1]*Q[1][1]},
                         {dR_*Q[0][0]*Q[1][0] + dT_*Q[0][1]*Q[1][1],
                          dR_*Q[1][0]*Q[1][0] + dT_*Q[1][1]*Q[1][1]}};
        return D;
    }

    Scalar massTransferCoefficient() const { return hM_; }
    Scalar airMoisture() const { return mAir_; }
    Scalar initialMoisture() const { return mInit_; }

    // unit radial direction at globalPos (from pith outward)
    GlobalPosition radialDirection(const GlobalPosition& globalPos) const
    {
        auto r = globalPos - pith_;
        r /= r.two_norm();
        return r;
    }

    // rotation Q whose columns are the local basis vectors expressed in the
    // global frame: Q = [e_R, e_T] with e_T the 90°-CCW rotation of e_R.
    // Q transforms vectors from the local (material) frame to the global frame.
    RotationMatrix rotation(const GlobalPosition& globalPos) const
    {
        const auto eR = radialDirection(globalPos);
        RotationMatrix Q;
        Q[0][0] =  eR[0]; Q[0][1] = -eR[1];
        Q[1][0] =  eR[1]; Q[1][1] =  eR[0];
        return Q;
    }

private:
    Scalar eR_, eT_, nuRT_, gRT_;
    Scalar alphaR_, alphaT_, mRef_;
    Scalar dR_, dT_, hM_, mAir_, mInit_;
    Scalar c11_, c22_, c12_, c66_;
    GlobalPosition pith_;
};

} // end namespace Dumux

#endif
