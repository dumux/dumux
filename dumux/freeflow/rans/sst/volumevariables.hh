// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::SSTMassVolumeVariables
 */
#ifndef DUMUX_RANS_SST_MASS_VOLUME_VARIABLES_HH
#define DUMUX_RANS_SST_MASS_VOLUME_VARIABLES_HH

#include <algorithm>
#include <cmath>

#include <dumux/common/parameters.hh>
#include <dumux/freeflow/turbulencemodel.hh>
#include <dumux/freeflow/navierstokes/mass/1p/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Volume variables for the single-phase Navier-Stokes mass balance fused with the
 *        two-equation SST (Menter 1994) turbulence closure for the turbulent kinetic energy k
 *        and the specific dissipation rate omega, fused as two extra equations on the mass
 *        sub-model rather than on a separate coupled sub-domain, following the same strategy
 *        already used for k-omega.
 *
 * Closure constants/functions ported near-verbatim from the deleted
 * releases/3.10:dumux/freeflow/rans/twoeq/sst/volumevariables.hh, with
 * stressTensorScalarProduct()/vorticityTensorScalarProduct()/wallDistance()/kinematicViscosity()
 * read from the problem (which forwards them to the momentum domain's
 * Dumux::RANSMomentumProblem), and the (lagged, see massproblem.hh) storedTurbulentKineticEnergy()/
 * storedDissipation()/storedCrossDiffusionGradientProduct() read from the mass problem's own
 * bookkeeping - the same pattern Dumux::KOmegaMassVolumeVariables already established. Unlike the
 * old code (which stored the two gradient vectors ∇k/∇ω separately), only their dot product is
 * stored, since that is the only quantity F1() and the cross-diffusion source term ever need -
 * the same simplification already made for k-omega's own cross-diffusion term.
 */
template<class Traits>
class SSTMassVolumeVariables
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

        const auto eIdx = scv.elementIndex();
        stressTensorScalarProduct_ = problem.stressTensorScalarProduct(eIdx);
        vorticityTensorScalarProduct_ = problem.vorticityTensorScalarProduct(eIdx);
        wallDistance_ = problem.wallDistance(eIdx);
        kinematicViscosity_ = problem.kinematicViscosity(eIdx);
        storedTurbulentKineticEnergy_ = problem.storedTurbulentKineticEnergy(eIdx);
        storedDissipation_ = problem.storedDissipation(eIdx);
        crossDiffusionGradientProduct_ = problem.storedCrossDiffusionGradientProduct(eIdx);

        // Cached here (rather than computed on demand from a Problem argument, the way
        // calculateEddyViscosity() itself needs one to branch on sstModelVersion()) so that a
        // plain no-argument dynamicEddyViscosity() getter can be exposed below - matching the
        // signature every other RANS model's own dynamicEddyViscosity() already has. This
        // matters beyond mere convention: dumux/multidomain/freeflow/couplingmanager_staggered.hh's
        // turbulentViscosity() (the sole path through which the momentum domain ever receives a
        // turbulent viscosity contribution) duck-types this exact zero-argument signature via
        // Dune::Std::is_detected; a Problem-argument-only overload is invisible to it and it
        // silently falls back to returning 0.0 - i.e. the momentum equation would run with
        // molecular viscosity only, no turbulence contribution at all. Recomputed here on every
        // update() call (i.e. every Newton iteration), so this is exactly as "live" as calling
        // calculateEddyViscosity(problem) directly would be - releases/3.10's optional
        // RANS.UseStoredEddyViscosity lagging opt-out is (deliberately) not ported, the same
        // simplification already made for the one-equation and k-omega models.
        dynamicEddyViscosity_ = calculateEddyViscosity(problem);

        // Nonisothermal only: see dumux/freeflow/rans/common/thermalconductivitymodel.hh's
        // class docs for why this must be done explicitly here rather than automatically.
        if constexpr (Traits::ModelTraits::enableEnergyBalance())
        {
            using FluidSystem = typename Traits::FluidSystem;
            static const Scalar turbulentPrandtlNumber = getParam<Scalar>("RANS.TurbulentPrandtlNumber", 1.0);
            this->lambdaEff_ = this->fluidThermalConductivity()
                + dynamicEddyViscosity_*FluidSystem::heatCapacity(this->fluidState(), 0)/turbulentPrandtlNumber;
        }
    }

    //! The turbulent kinetic energy primary variable.
    Scalar turbulentKineticEnergy() const
    { return this->priVar(Indices::turbulentKineticEnergyIdx); }

    //! The specific dissipation rate primary variable.
    Scalar dissipation() const
    { return this->priVar(Indices::dissipationIdx); }

    Scalar storedTurbulentKineticEnergy() const
    { return storedTurbulentKineticEnergy_; }

    Scalar storedDissipation() const
    { return storedDissipation_; }

    //! grad(k).grad(omega), from the problem's stored (lagged) gradients - needed by F1() and
    //! (only when positive) the cross-diffusion source term.
    Scalar storedCrossDiffusionGradientProduct() const
    { return crossDiffusionGradientProduct_; }

    Scalar stressTensorScalarProduct() const
    { return stressTensorScalarProduct_; }

    Scalar vorticityTensorScalarProduct() const
    { return vorticityTensorScalarProduct_; }

    Scalar wallDistance() const
    { return wallDistance_; }

    Scalar kinematicViscosity() const
    { return kinematicViscosity_; }

    //! \brief Returns the absolute value of the vorticity, Omega = sqrt(2*vorticityTensorScalarProduct).
    Scalar absoluteValueVorticity() const
    {
        using std::sqrt;
        return sqrt(2.0*vorticityTensorScalarProduct());
    }

    //! The (live, current-Newton-iterate) dynamic eddy viscosity, cached in update() - see its
    //! comment for why this needs to be a stored member with a no-argument getter rather than
    //! computed on demand from calculateEddyViscosity(problem) at every call site.
    Scalar dynamicEddyViscosity() const
    { return dynamicEddyViscosity_; }

    //! \brief Computes the dynamic eddy viscosity for the SST model, branching on the
    //!        (runtime-selected) SST model version. Only called from update() (see
    //!        dynamicEddyViscosity() above for the cached, no-argument accessor everything
    //!        else should use).
    template<class Problem>
    Scalar calculateEddyViscosity(const Problem& problem) const
    {
        using std::max;

        if (problem.sstModelVersion() == SSTModel::BSL)
            return turbulentKineticEnergy()/dissipation()*this->density();
        else
        {
            const Scalar dividend = a1SST()*turbulentKineticEnergy();
            const Scalar divisor = max(a1SST()*dissipation(), absoluteValueVorticity()*F2());
            return dividend/divisor*this->density();
        }
    }

    //! \brief Returns the transformation function F1 blending the near-wall (k-omega) and
    //!        far-field (transformed k-epsilon) closure-constant sets.
    //! \note Guarded against the transient storedDissipation()==0.0 that exists briefly between
    //! updateStaticWallProperties() (which only zero-initializes the lagged storage) and the
    //! first call to updateDynamicWallProperties() (which populates it from the real initial
    //! solution) - exercised once, when GridVariables::init() eagerly caches the t=0 volume
    //! variables in main.cc. Returning 0.0 there (a value F1 legitimately takes far from the
    //! wall) avoids a 0.0/0.0 NaN that would otherwise poison that cached state; it is never
    //! read by any solved time step (see the analogous guard for LowReKEpsilonMassProblem::yPlus()).
    Scalar F1() const
    {
        using std::max;
        using std::min;
        using std::sqrt;
        using std::tanh;

        if (storedDissipation() <= 0.0)
            return 0.0;

        Scalar positiveCrossDiffusion = 2.0*this->density()*sigmaOmega2()/storedDissipation()*storedCrossDiffusionGradientProduct();
        positiveCrossDiffusion = max(positiveCrossDiffusion, 1e-20);

        const Scalar possibleMin2 = (4.0*this->density()*sigmaOmega2()*storedTurbulentKineticEnergy())
                                  / (positiveCrossDiffusion*wallDistance()*wallDistance());

        const Scalar possibleMax1 = sqrt(storedTurbulentKineticEnergy())/(0.09*storedDissipation()*wallDistance());
        const Scalar possibleMax2 = (500.0*kinematicViscosity())/(wallDistance()*wallDistance()*storedDissipation());
        const Scalar possibleMin1 = max(possibleMax1, possibleMax2);

        const Scalar argument = min(possibleMin1, possibleMin2);
        return tanh(argument*argument*argument*argument);
    }

    //! \brief Returns the transformation function F2 used only by the eddy-viscosity limiter
    //!        (SST model version).
    Scalar F2() const
    {
        using std::max;
        using std::sqrt;
        using std::tanh;

        if (storedDissipation() <= 0.0)
            return 0.0;

        const Scalar possibleMax1 = 2.0*sqrt(storedTurbulentKineticEnergy())/(0.09*storedDissipation()*wallDistance());
        const Scalar possibleMax2 = (500.0*kinematicViscosity())/(wallDistance()*wallDistance()*storedDissipation());
        const Scalar argument = max(possibleMax1, possibleMax2);
        return tanh(argument*argument);
    }

    //! \name Closure constants (near-wall "1" set, far-field "2" set, and the two blended sets)
    // \{
    Scalar betaOmega() const { return 0.0708; }

    Scalar sigmaK1BSL() const { return 0.5; }
    Scalar sigmaK1SST() const { return 0.85; }
    Scalar sigmaK2() const { return 1.0; }

    Scalar sigmaOmega1BSL() const { return 0.5; }
    Scalar sigmaOmega1SST() const { return 0.5; }
    Scalar sigmaOmega2() const { return 0.856; }

    Scalar beta1BSL() const { return 0.0750; }
    Scalar beta1SST() const { return 0.0750; }
    Scalar beta2() const { return 0.0820; }

    Scalar betaStar1BSL() const { return 0.09; }
    Scalar betaStar1SST() const { return 0.09; }
    Scalar betaStar2() const { return 0.09; }

    Scalar kappa1BSL() const { return 0.41; }
    Scalar kappa1SST() const { return 0.41; }
    Scalar kappa2() const { return 0.41; }

    Scalar gamma1BSL() const
    {
        using std::sqrt;
        return beta1BSL()/betaStar1BSL() - sigmaOmega1BSL()*kappa1BSL()*kappa1BSL()/sqrt(betaStar1BSL());
    }
    Scalar gamma1SST() const
    {
        using std::sqrt;
        return beta1SST()/betaStar1SST() - sigmaOmega1SST()*kappa1SST()*kappa1SST()/sqrt(betaStar1SST());
    }
    Scalar gamma2() const
    {
        using std::sqrt;
        return beta2()/betaStar2() - sigmaOmega2()*kappa2()*kappa2()/sqrt(betaStar2());
    }

    Scalar a1SST() const { return 0.31; }

    Scalar sigmaKSST() const { return F1()*sigmaK1SST() + (1.0 - F1())*sigmaK2(); }
    Scalar sigmaOmegaSST() const { return F1()*sigmaOmega1SST() + (1.0 - F1())*sigmaOmega2(); }
    Scalar betaSST() const { return F1()*beta1SST() + (1.0 - F1())*beta2(); }
    Scalar betaStarSST() const { return F1()*betaStar1SST() + (1.0 - F1())*betaStar2(); }
    Scalar gammaSST() const { return F1()*gamma1SST() + (1.0 - F1())*gamma2(); }

    Scalar sigmaKBSL() const { return F1()*sigmaK1BSL() + (1.0 - F1())*sigmaK2(); }
    Scalar sigmaOmegaBSL() const { return F1()*sigmaOmega1BSL() + (1.0 - F1())*sigmaOmega2(); }
    Scalar betaBSL() const { return F1()*beta1BSL() + (1.0 - F1())*beta2(); }
    Scalar betaStarBSL() const { return F1()*betaStar1BSL() + (1.0 - F1())*betaStar2(); }
    Scalar gammaBSL() const { return F1()*gamma1BSL() + (1.0 - F1())*gamma2(); }
    // \}

    //! \name Model-version-dependent accessors, used by the diffusive flux (localresidual.hh)
    //!       so it does not need to branch on problem.sstModelVersion() itself.
    // \{
    template<class Problem>
    Scalar sigmaK(const Problem& problem) const
    { return problem.sstModelVersion() == SSTModel::SST ? sigmaKSST() : sigmaKBSL(); }

    template<class Problem>
    Scalar sigmaOmega(const Problem& problem) const
    { return problem.sstModelVersion() == SSTModel::SST ? sigmaOmegaSST() : sigmaOmegaBSL(); }
    // \}

protected:
    Scalar stressTensorScalarProduct_ = 0.0;
    Scalar vorticityTensorScalarProduct_ = 0.0;
    Scalar wallDistance_ = 0.0;
    Scalar kinematicViscosity_ = 0.0;
    Scalar storedTurbulentKineticEnergy_ = 0.0;
    Scalar storedDissipation_ = 0.0;
    Scalar crossDiffusionGradientProduct_ = 0.0;
    Scalar dynamicEddyViscosity_ = 0.0;
};

} // end namespace Dumux

#endif
