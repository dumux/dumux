// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::OneEqMassVolumeVariables
 */
#ifndef DUMUX_RANS_ONEEQ_MASS_VOLUME_VARIABLES_HH
#define DUMUX_RANS_ONEEQ_MASS_VOLUME_VARIABLES_HH

#include <algorithm>
#include <cmath>

#include <dumux/freeflow/navierstokes/mass/1p/volumevariables.hh>
#include <dumux/freeflow/rans/common/thermalconductivitymodel.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Volume variables for the single-phase Navier-Stokes mass balance fused with the
 *        one-equation (Spalart-Allmaras) turbulence transport equation for the working
 *        viscosity ν̃ (see turbulenceequations.md \S7 for the physics, and
 *        whatisimplemented.md/proposedimplementation.md for why ν̃ lives here - as a second
 *        equation on the mass sub-model - rather than on a separate coupled sub-domain).
 *
 * The Spalart-Allmaras closure functions here are ported near-verbatim from the deleted
 * releases/3.10:dumux/freeflow/rans/oneeq/volumevariables.hh, with wallDistance()/
 * vorticityTensorScalarProduct() and the (lagged, see problem.hh) storedViscosityTildeGradient()
 * read from the problem (which forwards the former two to the momentum domain's
 * Dumux::RANSMomentumProblem, and computes the latter itself, see
 * dumux/freeflow/rans/oneeq/massproblem.hh).
 */
template<class Traits>
class OneEqMassVolumeVariables
: public NavierStokesMassOnePVolumeVariables<Traits>
{
    using ParentType = NavierStokesMassOnePVolumeVariables<Traits>;
    using Scalar = typename Traits::PrimaryVariables::value_type;

public:
    using Indices = typename Traits::ModelTraits::Indices;

    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        wallDistance_ = problem.wallDistance(scv.elementIndex());
        vorticityTensorScalarProduct_ = problem.vorticityTensorScalarProduct(scv.elementIndex());
        viscosityTildeGradientSquared_ = problem.storedViscosityTildeGradient(scv.elementIndex()).two_norm2();

        // Nonisothermal only: see dumux/freeflow/rans/common/thermalconductivitymodel.hh's
        // class docs for why this must be done explicitly here rather than automatically.
        if constexpr (Traits::ModelTraits::enableEnergyBalance())
            this->lambdaEff_ = RANSThermalConductivityModel::turbulentEffectiveThermalConductivity(*this);
    }

    //! The working (Spalart-Allmaras) viscosity primary variable.
    Scalar viscosityTilde() const
    { return this->priVar(Indices::viscosityTildeIdx); }

    Scalar kinematicViscosity() const
    { return this->viscosity()/this->density(); }

    Scalar wallDistance() const
    { return wallDistance_; }

    Scalar karmanConstant() const
    { return 0.41; }

    //! The dynamic eddy viscosity mu_t = rho * nu~ * fv1 - the quantity the momentum domain
    //! reads through the coupling manager's turbulentViscosity()/effectiveViscosity().
    Scalar dynamicEddyViscosity() const
    { return viscosityTilde()*fv1()*this->density(); }

    Scalar viscosityRatio() const
    { return viscosityTilde()/kinematicViscosity(); }

    Scalar fv1() const
    {
        const Scalar chi3 = viscosityRatio()*viscosityRatio()*viscosityRatio();
        return chi3/(chi3 + cv1()*cv1()*cv1());
    }

    Scalar fv2() const
    { return 1.0 - viscosityRatio()/(1.0 + viscosityRatio()*fv1()); }

    Scalar fw() const
    {
        using std::pow;
        const Scalar g_ = g();
        const Scalar cw3_6 = cw3()*cw3()*cw3()*cw3()*cw3()*cw3();
        const Scalar g_6 = g_*g_*g_*g_*g_*g_;
        return g_*pow((1.0 + cw3_6)/(g_6 + cw3_6), 1.0/6.0);
    }

    Scalar g() const
    {
        const Scalar r_ = r();
        const Scalar r_6 = r_*r_*r_*r_*r_*r_;
        return r_ + cw2()*(r_6 - r_);
    }

    Scalar r() const
    {
        using std::min;
        using std::max;
        // Guard against stressTensorScalarProductTilde() == 0 exactly (e.g. at a channel
        // centerline, where the mean velocity profile's symmetry makes the vorticity - and
        // hence, in the sBar<-c2*omega branch, the whole modified strain rate - exactly zero):
        // without this floor, r() -> infinity, cascading through g()/fw() into an Inf*0 = NaN.
        const Scalar sTilde = max(stressTensorScalarProductTilde(), 1e-10);
        return min(10.0, viscosityTilde()/(sTilde*karmanConstant()*karmanConstant()*wallDistance()*wallDistance()));
    }

    //! Modified mean effective strain rate (Allmaras/Johnson/Spalart limiter form).
    Scalar stressTensorScalarProductTilde() const
    {
        const Scalar sBar = viscosityTilde()*fv2()/(karmanConstant()*karmanConstant()*wallDistance()*wallDistance());
        const Scalar omega = vorticityMagnitude();

        if (sBar < -c2()*omega)
            return omega + (omega*(c2()*c2()*omega + c3()*sBar))/((c3() - 2.0*c2())*omega - sBar);
        else
            return omega + sBar;
    }

    Scalar vorticityMagnitude() const
    {
        using std::sqrt;
        return sqrt(2.0*vorticityTensorScalarProduct_);
    }

    //! |grad(nu~)|^2 at this element, from the problem's stored (lagged) gradient, needed for
    //! the cross-diffusion source term (cb2/sigma)*rho*|grad(nu~)|^2.
    Scalar viscosityTildeGradientSquared() const
    { return viscosityTildeGradientSquared_; }

    Scalar c2() const { return 0.7; }
    Scalar c3() const { return 0.9; }
    Scalar sigma() const { return 2.0/3.0; }
    Scalar cb1() const { return 0.1355; }
    Scalar cb2() const { return 0.622; }
    Scalar cv1() const { return 7.1; }
    Scalar cw2() const { return 0.3; }
    Scalar cw3() const { return 2.0; }

    Scalar cw1() const
    { return cb1()/(karmanConstant()*karmanConstant()) + (1.0 + cb2())/sigma(); }

protected:
    Scalar wallDistance_ = 0.0;
    Scalar vorticityTensorScalarProduct_ = 0.0;
    Scalar viscosityTildeGradientSquared_ = 0.0;
};

} // end namespace Dumux

#endif
