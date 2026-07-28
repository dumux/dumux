// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::KOmegaMassVolumeVariables
 */
#ifndef DUMUX_RANS_KOMEGA_MASS_VOLUME_VARIABLES_HH
#define DUMUX_RANS_KOMEGA_MASS_VOLUME_VARIABLES_HH

#include <algorithm>
#include <cmath>

#include <dumux/freeflow/navierstokes/mass/1p/volumevariables.hh>
#include <dumux/freeflow/rans/common/thermalconductivitymodel.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Volume variables for the single-phase Navier-Stokes mass balance fused with the
 *        two-equation k-omega (Wilcox 2008) turbulence closure for the turbulent kinetic
 *        energy k and the specific dissipation rate omega (see turbulenceequations.md \S8 for
 *        the physics, and whatisimplemented.md/proposedimplementation.md for why k/omega live
 *        here - as two extra equations on the mass sub-model - rather than on a separate
 *        coupled sub-domain, following the same strategy validated for the one-equation model).
 *
 * Closure constants/functions ported near-verbatim from the deleted
 * releases/3.10:dumux/freeflow/rans/twoeq/komega/volumevariables.hh, with
 * stressTensorScalarProduct() read from the problem (which forwards it to the momentum
 * domain's Dumux::RANSMomentumProblem), and the (lagged, see massproblem.hh)
 * storedTurbulentKineticEnergyGradient()/storedDissipationGradient() read from the mass
 * problem's own bookkeeping.
 */
template<class Traits>
class KOmegaMassVolumeVariables
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

        stressTensorScalarProduct_ = problem.stressTensorScalarProduct(scv.elementIndex());
        crossDiffusionGradientProduct_ = problem.storedCrossDiffusionGradientProduct(scv.elementIndex());

        // Nonisothermal only: overwrite the molecular-only lambdaEff_ that
        // ParentType::update() already set (via a call site that cannot see
        // dynamicEddyViscosity(), see dumux/freeflow/rans/common/thermalconductivitymodel.hh's
        // class docs) with the correct eddy-augmented value, now that dynamicEddyViscosity() -
        // a pure function of the primary variables ParentType::update() just set up - is
        // reachable through *this.
        if constexpr (Traits::ModelTraits::enableEnergyBalance())
            this->lambdaEff_ = RANSThermalConductivityModel::turbulentEffectiveThermalConductivity(*this);
    }

    //! The turbulent kinetic energy primary variable.
    Scalar turbulentKineticEnergy() const
    { return this->priVar(Indices::turbulentKineticEnergyIdx); }

    //! The specific dissipation rate primary variable.
    Scalar dissipation() const
    { return this->priVar(Indices::dissipationIdx); }

    Scalar stressTensorScalarProduct() const
    { return stressTensorScalarProduct_; }

    //! grad(k).grad(omega), from the problem's stored (lagged) gradients, needed for the
    //! cross-diffusion source term (only used - by the caller - when this is positive).
    Scalar storedCrossDiffusionGradientProduct() const
    { return crossDiffusionGradientProduct_; }

    //! The dynamic eddy viscosity mu_t = rho*k/max(omega, limiter) - the quantity the momentum
    //! domain reads through the coupling manager's turbulentViscosity()/effectiveViscosity().
    //! The Wilcox2008 dissipation limiter is always enabled (releases/3.10's
    //! RANS.UseStoredEddyViscosity-style runtime opt-out is not ported, matching the
    //! simplification already made for the one-equation model, see whatisimplemented.md).
    Scalar dynamicEddyViscosity() const
    {
        using std::sqrt;
        using std::max;

        const Scalar limitedDissipation = (7.0/8.0)*sqrt(2.0*stressTensorScalarProduct()/betaK());
        return turbulentKineticEnergy()/max(dissipation(), limitedDissipation)*this->density();
    }

    Scalar alpha() const { return 0.52; }
    Scalar sigmaK() const { return 0.6; }
    Scalar sigmaOmega() const { return 0.5; }
    Scalar betaK() const { return 0.09; }
    Scalar betaOmega() const { return 0.0708; }

protected:
    Scalar stressTensorScalarProduct_ = 0.0;
    Scalar crossDiffusionGradientProduct_ = 0.0;
};

} // end namespace Dumux

#endif
