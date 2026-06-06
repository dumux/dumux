// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Spatial parameters for the Cryer sphere consolidation benchmark.
 *
 * Material: Neo-Hookean hyperelastic porous medium with:
 * - Shear modulus G = 1.5 MPa
 * - Bulk modulus  K = 1.0 MPa  (Poisson's ratio nu = 0; W_vol uses K directly)
 * - Biot coefficient alpha_B = 1
 * - Storage coefficient Sp = 0 (incompressible solid/fluid)
 * - Initial porosity phi_0 = 0.5
 * - Constant or Kozeny–Carman permeability
 */
#ifndef DUMUX_CRYER_SPATIALPARAMS_HH
#define DUMUX_CRYER_SPATIALPARAMS_HH

#include <cmath>
#include <string>
#include <dumux/common/parameters.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>

namespace Dumux {

/*!
 * \brief Spatial parameters for the three subdomains of the Cryer benchmark.
 *
 * All three subdomains share the same spatial params object (the same underlying
 * material parameters apply everywhere).
 */
template<class GridGeometry, class Scalar>
class CryerSpatialParams
{
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

public:
    CryerSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_(gridGeometry)
    {
        G_     = getParam<Scalar>("SpatialParams.ShearModulus");
        K_     = getParam<Scalar>("SpatialParams.BulkModulus");
        alphaB_ = getParam<Scalar>("SpatialParams.BiotCoefficient", 1.0);
        Sp_     = getParam<Scalar>("SpatialParams.StorageCoefficient", 0.0);
        phi0_   = getParam<Scalar>("SpatialParams.InitialPorosity");
        k0_     = getParam<Scalar>("SpatialParams.InitialPermeability");
        muf_    = getParam<Scalar>("SpatialParams.FluidViscosity");

        const std::string permModel = getParam<std::string>("SpatialParams.PermeabilityModel", "Constant");
        useKozenyCarman_ = (permModel == "KozenyCarman");

        // Kozeny-Carman particle diameter: dp = sqrt(180*k0*(1-phi0)^2/phi0^3)
        dp_ = getParam<Scalar>("SpatialParams.ParticleDiameter",
            std::sqrt(180.0 * k0_ * (1.0-phi0_)*(1.0-phi0_) / (phi0_*phi0_*phi0_)));
    }

    //! Shear modulus G [Pa]
    Scalar shearModulus() const { return G_; }

    //! Bulk modulus K [Pa]
    Scalar bulkModulus() const { return K_; }

    //! Biot coefficient alpha_B [-]
    template<class SubControlVolume>
    Scalar biotCoefficient(const Element&, const SubControlVolume&) const { return alphaB_; }

    //! Storage coefficient Sp [1/Pa]
    template<class SubControlVolume>
    Scalar storageCoefficient(const Element&, const SubControlVolume&) const { return Sp_; }

    //! Initial porosity phi_0 [-]
    Scalar initialPorosity() const { return phi0_; }

    //! Fluid dynamic viscosity mu_f [Pa·s]
    template<class FVElementGeometry>
    Scalar fluidViscosity(const Element&, const FVElementGeometry&) const { return muf_; }

    /*!
     * \brief Permeability at a given point as a function of the Jacobian determinant J.
     *
     * If the Kozeny–Carman model is selected:
     *   phi(J) = 1 - (1-phi0)/J
     *   k(J)   = dp^2/180 * phi(J)^3 / (1-phi(J))^2
     *
     * Otherwise: k = k0 (constant).
     */
    template<class FVElementGeometry>
    Scalar permeabilityAtPoint(const Element&, const FVElementGeometry&,
                               const typename GridGeometry::GlobalCoordinate&,
                               Scalar J) const
    {
        if (!useKozenyCarman_)
            return k0_;

        // updated porosity from solid volume conservation
        const Scalar phi = 1.0 - (1.0 - phi0_) / J;
        // Kozeny–Carman
        return dp_*dp_ / 180.0 * phi*phi*phi / ((1.0-phi)*(1.0-phi));
    }

    const GridGeometry& gridGeometry() const { return *gridGeometry_; }

private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
    Scalar G_, K_, alphaB_, Sp_, phi0_, k0_, muf_, dp_;
    bool useKozenyCarman_ = false;
};

} // end namespace Dumux
#endif
