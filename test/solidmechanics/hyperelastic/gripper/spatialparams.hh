// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Spatial parameters for the gripper hydrogel actuator.
 *
 * The domain is a bilayer cantilever beam:
 * - **Active layer** (top half, y > H/2): soft hydrogel with low stiffness and
 *   high permeability; osmotic swelling is driven here.
 * - **Passive layer** (bottom half, y ≤ H/2): stiff passive material with low
 *   permeability; provides the bending resistance.
 *
 * Permeability follows a Kozeny-Carman-inspired exponential law:
 * \f[ K(J) = K_0 \exp\!\left(M \frac{J-1}{J}\right) \f]
 *
 * Osmotic pressure from Flory-Huggins theory (used by the pressure problem for BCs):
 * \f[ \Pi(J) = -\frac{RT}{V_m}\left[\ln(1-\phi) + \phi + \chi\phi^2\right],
 *     \quad \phi = \phi_0/J \f]
 */
#ifndef DUMUX_GRIPPER_SPATIAL_PARAMS_HH
#define DUMUX_GRIPPER_SPATIAL_PARAMS_HH

#include <cmath>

#include <dumux/common/fvspatialparams.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \brief Spatial parameters shared by both gripper subdomains.
 *
 * \tparam GridGeometry The grid-geometry type of either subdomain (both share the same grid).
 * \tparam Scalar       The scalar type.
 */
template<class GridGeometry, class Scalar>
class GripperSpatialParams
: public FVSpatialParams<GridGeometry, Scalar, GripperSpatialParams<GridGeometry, Scalar>>
{
    using ThisType   = GripperSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView   = typename GridGeometry::GridView;
    using Element    = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    GripperSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // Bilayer geometry: active layer is y > H/2.
        yMax_ = this->gridGeometry().bBoxMax()[1];
        yMin_ = this->gridGeometry().bBoxMin()[1];
        layerInterface_ = 0.5 * (yMin_ + yMax_);
        xMax_ = this->gridGeometry().bBoxMax()[0];
        xMin_ = this->gridGeometry().bBoxMin()[0];
        activeRegionStart_ = getParam<Scalar>("Problem.ActiveRegionStart", Scalar(0.2));

        // Stiffness (Lamé parameters) for each layer.
        muActive_ = getParam<Scalar>("SpatialParams.ShearModulusActive");
        lambdaActive_ = getParam<Scalar>("SpatialParams.FirstLameActive");
        muPassive_ = getParam<Scalar>("SpatialParams.ShearModulusPassive");
        lambdaPassive_ = getParam<Scalar>("SpatialParams.FirstLamePassive");

        // Permeability.
        K0active_ = getParam<Scalar>("SpatialParams.PermeabilityActive");
        K0passive_ = getParam<Scalar>("SpatialParams.PermeabilityPassive");
        permExpM_ = getParam<Scalar>("SpatialParams.PermeabilityExponent", Scalar(1.0));

        // Fluid properties.
        muFluid_ = getParam<Scalar>("SpatialParams.FluidViscosity");

        // Flory-Huggins osmotic pressure parameters.
        phi0_ = getParam<Scalar>("SpatialParams.SolidVolumeFraction");
        chi_ = getParam<Scalar>("SpatialParams.FloryHugginsParameter");
        R_ = Scalar(8.314);   // J/(mol·K)
        T_ = getParam<Scalar>("SpatialParams.Temperature", Scalar(298.15));
        Vm_ = getParam<Scalar>("SpatialParams.MolarVolume");
    }

    //! True if the element belongs to the active (swelling) layer.
    bool isActiveLayer(const Element& element) const
    {
        const auto c = element.geometry().center();
        return c[1] > layerInterface_
               && c[0] - xMin_ > activeRegionStart_ * (xMax_ - xMin_);
    }

    //! Neo-Hookean shear modulus \f$\mu\f$ for the element's layer.
    Scalar shearModulus(const Element& element) const
    { return isActiveLayer(element) ? muActive_ : muPassive_; }

    //! Neo-Hookean first Lamé parameter \f$\lambda\f$ for the element's layer.
    Scalar firstLame(const Element& element) const
    { return isActiveLayer(element) ? lambdaActive_ : lambdaPassive_; }

    /*!
     * \brief Permeability \f$K(J)\f$ for the element's layer.
     *
     * \f[ K(J) = K_0 \exp\!\left(M \frac{J-1}{J}\right) \f]
     */
    Scalar permeability(const Element& element, Scalar J) const
    {
        const Scalar K0 = isActiveLayer(element) ? K0active_ : K0passive_;
        return K0 * std::exp(permExpM_ * (J - 1.0) / std::max(J, Scalar(1e-8)));
    }

    //! Fluid (solvent) dynamic viscosity \f$\mu_f\f$.
    Scalar fluidViscosity() const { return muFluid_; }

    //! Reference solid volume fraction \f$\phi_0\f$ (used for osmotic BC).
    Scalar solidVolumeFraction() const { return phi0_; }

    /*!
     * \brief Flory-Huggins osmotic pressure \f$\Pi(J)\f$ for the active layer.
     *
     * \f[
     *   \Pi(J) = -\frac{RT}{V_m}\left[\ln(1-\phi) + \phi + \chi\phi^2\right],
     *   \quad \phi = \phi_0 / J
     * \f]
     *
     * Returned as a positive value when the hydrogel wants to absorb solvent.
     */
    Scalar osmoticPressure(Scalar J) const
    {
        const Scalar phi = phi0_ / std::max(J, Scalar(1e-8));
        const Scalar phi_clip = std::min(phi, Scalar(1.0 - 1e-8));
        return -(R_ * T_ / Vm_) * (std::log(1.0 - phi_clip) + phi_clip + chi_ * phi_clip * phi_clip);
    }

    /*!
     * \brief Convenience: compute Lamé parameters from Young's modulus and Poisson's ratio.
     * \return {mu, lambda} pair.
     */
    static std::pair<Scalar, Scalar> fromYoungsPoisson(Scalar E, Scalar nu)
    {
        const Scalar mu     = E / (2.0 * (1.0 + nu));
        const Scalar lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        return {mu, lambda};
    }

private:
    Scalar yMin_, yMax_, layerInterface_, xMin_, xMax_, activeRegionStart_;
    Scalar muActive_, lambdaActive_;
    Scalar muPassive_, lambdaPassive_;
    Scalar K0active_, K0passive_, permExpM_;
    Scalar muFluid_;
    Scalar phi0_, chi_, R_, T_, Vm_;
};

} // end namespace Dumux

#endif // DUMUX_GRIPPER_SPATIAL_PARAMS_HH
