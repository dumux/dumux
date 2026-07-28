// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::LowReKEpsilonMassVolumeVariables
 */
#ifndef DUMUX_RANS_LOWREKEPSILON_MASS_VOLUME_VARIABLES_HH
#define DUMUX_RANS_LOWREKEPSILON_MASS_VOLUME_VARIABLES_HH

#include <cmath>

#include <dumux/freeflow/navierstokes/mass/1p/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Volume variables for the single-phase Navier-Stokes mass balance fused with the
 *        two-equation, low-Reynolds-number k-epsilon (Chien 1982) turbulence closure for the
 *        turbulent kinetic energy k and the transported dissipation epsilon-tilde (see
 *        whatisimplemented.md/proposedimplementation.md for why k/epsilon-tilde live here - as
 *        two extra equations on the mass sub-model - rather than on a separate coupled
 *        sub-domain, following the same strategy validated for the one-equation and k-omega
 *        models).
 *
 * Unlike the wall-function k-epsilon model (dumux/freeflow/rans/kepsilon/), this model damps the
 * eddy viscosity and the epsilon-tilde production/destruction terms analytically as functions of
 * y+ and the turbulent Reynolds number Re_t, so that k=0, epsilon-tilde=0 ordinary Dirichlet
 * conditions at the wall are enough - no near-wall-region/matching-point bookkeeping is needed at
 * all. Closure constants/functions ported near-verbatim from the deleted
 * releases/3.10:dumux/freeflow/rans/twoeq/lowrekepsilon/volumevariables.hh, with
 * stressTensorScalarProduct()/yPlus() read from the problem (which forwards them to the momentum
 * domain's Dumux::RANSMomentumProblem).
 */
template<class Traits>
class LowReKEpsilonMassVolumeVariables
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
        yPlus_ = problem.yPlus(scv.elementIndex());
        wallDistance_ = problem.wallDistance(scv.elementIndex());
    }

    //! The turbulent kinetic energy primary variable.
    Scalar turbulentKineticEnergy() const
    { return this->priVar(Indices::turbulentKineticEnergyIdx); }

    //! The transported dissipation (epsilon-tilde) primary variable.
    Scalar dissipationTilde() const
    { return this->priVar(Indices::dissipationIdx); }

    Scalar stressTensorScalarProduct() const
    { return stressTensorScalarProduct_; }

    //! The wall unit y+, from the problem's (lagged) local wall-shear-velocity bookkeeping.
    Scalar yPlus() const
    { return yPlus_; }

    //! The turbulent Reynolds number Re_t = k^2/(nu*epsilonTilde).
    Scalar reT() const
    { return turbulentKineticEnergy()*turbulentKineticEnergy()/(this->kinematicViscosity()*dissipationTilde()); }

    Scalar kinematicViscosity() const
    { return this->viscosity()/this->density(); }

    Scalar cMu() const { return 0.09; }
    Scalar sigmaK() const { return 1.0; }
    Scalar sigmaEpsilon() const { return 1.3; }
    Scalar cOneEpsilon() const { return 1.35; }
    Scalar cTwoEpsilon() const { return 1.8; }

    //! Chien's damping function f_mu, damping the eddy viscosity to zero at the wall.
    Scalar fMu() const
    {
        using std::exp;
        return 1.0 - exp(-0.0115*yPlus());
    }

    //! Chien's f_1 (no additional damping of the production term in this variant).
    Scalar fOne() const
    { return 1.0; }

    //! Chien's f_2, damping the destruction term where Re_t is small (near the wall).
    Scalar fTwo() const
    {
        using std::exp;
        const Scalar ret = reT();
        return 1.0 - 0.22*exp(-(ret*ret)/36.0);
    }

    //! The dynamic eddy viscosity, damped by f_mu: mu_t = rho*cMu*fMu*k^2/epsilonTilde.
    Scalar dynamicEddyViscosity() const
    { return cMu()*fMu()*turbulentKineticEnergy()*turbulentKineticEnergy()/dissipationTilde()*this->density(); }

    //! The extra sink D = 2*nu*k/y^2 in the k-equation, letting epsilon-tilde asymptote to the
    //! true dissipation epsilon = epsilonTilde + D at the wall.
    Scalar dValue() const
    {
        const Scalar y = wallDistance_;
        return 2.0*kinematicViscosity()*turbulentKineticEnergy()/(y*y);
    }

    //! The extra source/sink E in the epsilon-tilde equation (see dValue()).
    Scalar eValue() const
    {
        using std::exp;
        const Scalar y = wallDistance_;
        return -2.0*kinematicViscosity()*dissipationTilde()/(y*y)*exp(-0.5*yPlus());
    }

protected:
    Scalar stressTensorScalarProduct_ = 0.0;
    Scalar yPlus_ = 0.0;
    Scalar wallDistance_ = 0.0;
};

} // end namespace Dumux

#endif
