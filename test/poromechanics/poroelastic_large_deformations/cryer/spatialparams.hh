// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Spatial parameters for the Cryer sphere consolidation benchmark.
 */
#ifndef DUMUX_TEST_POROELASTIC_LARGE_DEF_CRYER_SPATIALPARAMS_HH
#define DUMUX_TEST_POROELASTIC_LARGE_DEF_CRYER_SPATIALPARAMS_HH

#include <cmath>
#include <memory>
#include <string>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief Spatial parameters for the three subdomains of the Cryer benchmark.
 */
template<class GridDiscretization, class Scalar>
class CryerSpatialParams
{
    using GlobalPosition = typename GridDiscretization::GlobalCoordinate;
    using Element = typename GridDiscretization::GridView::template Codim<0>::Entity;

public:
    CryerSpatialParams(std::shared_ptr<const GridDiscretization> gridDiscretization)
    : gridDiscretization_(gridDiscretization)
    {
        G_ = getParam<Scalar>("SpatialParams.ShearModulus");
        K_ = getParam<Scalar>("SpatialParams.BulkModulus");
        alphaB_ = getParam<Scalar>("SpatialParams.BiotCoefficient", 1.0);
        Sp_ = getParam<Scalar>("SpatialParams.StorageCoefficient", 0.0);
        phi0_ = getParam<Scalar>("SpatialParams.InitialPorosity");
        k0_ = getParam<Scalar>("SpatialParams.InitialPermeability");
        muf_ = getParam<Scalar>("SpatialParams.FluidViscosity");

        const std::string permModel = getParam<std::string>("SpatialParams.PermeabilityModel", "Constant");
        useKozenyCarman_ = (permModel == "KozenyCarman");

        // such that the Kozeny-Carman permeability at J = 1 is k0
        dp_ = getParam<Scalar>("SpatialParams.ParticleDiameter",
            std::sqrt(180.0 * k0_ * (1.0-phi0_)*(1.0-phi0_) / (phi0_*phi0_*phi0_)));
    }

    //! Shear modulus G [Pa]
    Scalar shearModulus() const { return G_; }

    //! Bulk modulus K [Pa]
    Scalar bulkModulus() const { return K_; }

    //! Biot coefficient alpha_B [-]
    template<class ElementDiscretization>
    Scalar biotCoefficient(const Element&, const ElementDiscretization&, const GlobalPosition&) const
    { return alphaB_; }

    //! Storage coefficient Sp = 1/M [1/Pa]
    template<class ElementDiscretization>
    Scalar storageCoefficient(const Element&, const ElementDiscretization&, const GlobalPosition&) const
    { return Sp_; }

    //! Initial porosity phi_0 [-]
    Scalar initialPorosity() const { return phi0_; }

    //! Fluid dynamic viscosity mu_f [Pa·s]
    template<class ElementDiscretization>
    Scalar fluidViscosity(const Element&, const ElementDiscretization&) const { return muf_; }

    /*!
     * \brief Constant or Kozeny–Carman permeability, the latter with the porosity
     *        \f$ \phi = 1 - (1 - \phi_0)/J \f$ from the conservation of the solid volume
     */
    template<class ElementDiscretization>
    Scalar permeability(const Element&, const ElementDiscretization&,
                               const typename GridDiscretization::GlobalCoordinate&,
                               Scalar J) const
    {
        if (!useKozenyCarman_)
            return k0_;

        const Scalar phi = 1.0 - (1.0 - phi0_) / J;
        return dp_*dp_ / 180.0 * phi*phi*phi / ((1.0-phi)*(1.0-phi));
    }

    const GridDiscretization& gridDiscretization() const { return *gridDiscretization_; }

private:
    std::shared_ptr<const GridDiscretization> gridDiscretization_;
    Scalar G_, K_, alphaB_, Sp_, phi0_, k0_, muf_, dp_;
    bool useKozenyCarman_ = false;
};

} // end namespace Dumux
#endif
