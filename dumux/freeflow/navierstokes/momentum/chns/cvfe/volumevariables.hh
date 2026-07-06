// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesMomentumCHNSCVFEVolumeVariables
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_CHNS_CVFE_VOLUME_VARIABLES_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_CHNS_CVFE_VOLUME_VARIABLES_HH

#include <dune/common/fvector.hh>
#include <dumux/common/concepts/ipdata_.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Volume variables for the co-located momentum + Cahn-Hilliard CVFE model.
 *
 * The primary variables are [ v_0 ... v_{dim-1}, c, mu ] (numEq = dim+2). Unlike the plain
 * momentum volume variables (where velocity() returns the whole priVars vector), here
 * velocity() returns ONLY the first dim components (a GlobalPosition), because the momentum
 * flux helpers static_assert NumEqVector::dimension == dimWorld and interpolate velocity as a
 * dim-vector. The phase field c and chemical potential mu are the last two components.
 *
 * Density and viscosity are MIXTURE quantities rho(c), nu(c) evaluated from the LOCAL phase
 * field via the problem's material laws (problem.mixtureDensity/mixtureViscosity). They are
 * exposed through density()/effectiveViscosity() so the (forked) momentum flux/residual read
 * them from the volvars, which is required for a correct monolithic Jacobian: the numeric
 * Jacobian deflects c per dof and rebuilds the volvars, so a volvars-local rho(c) responds to
 * the deflection whereas a coupling-manager lookup of a stored solution would not.
 */
template <class Traits>
class NavierStokesMomentumCHNSCVFEVolumeVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using ModelTraits = typename Traits::ModelTraits;

    static constexpr int dim = ModelTraits::dim();
    static_assert(Traits::PrimaryVariables::dimension == dim + 2,
                  "CHNS momentum model expects numEq = dim + 2 (velocity, c, mu)");

public:
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using Indices = typename ModelTraits::Indices;
    using GlobalPosition = Dune::FieldVector<Scalar, dim>;

    /*!
     * \brief Update all quantities for a given local dof (new-variables interface).
     * \note The VariablesAdapter calls this per local dof with the LocalDofIpData, so the priVars
     *       are read via ipData.localDofIndex() (covers vertex AND P2 edge dofs).
     */
    template<class ElementSolution, class Problem, class FVElementGeometry, Concept::LocalDofIpData IpData>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const FVElementGeometry& fvGeometry,
                const IpData& ipData)
    {
        priVars_ = elemSol[ipData.localDofIndex()];
        extrusionFactor_ = problem.spatialParams().extrusionFactor(fvGeometry, ipData, elemSol);

        // mixture material from the LOCAL phase field (needed for a correct monolithic Jacobian)
        density_ = problem.mixtureDensity(phaseField());
        viscosity_ = problem.mixtureViscosity(phaseField());
    }

    //! velocity is the first dim components of the primary variables
    GlobalPosition velocity() const
    {
        GlobalPosition v(0.0);
        for (int i = 0; i < dim; ++i)
            v[i] = priVars_[Indices::velocity(i)];
        return v;
    }

    Scalar velocity(const int dirIdx) const
    { return priVars_[Indices::velocity(dirIdx)]; }

    //! the phase field c (= +1 inside the bubble, -1 outside; test convention)
    Scalar phaseField() const
    { return priVars_[Indices::phaseFieldIdx]; }

    //! the chemical potential mu
    Scalar chemicalPotential() const
    { return priVars_[Indices::chemicalPotentialIdx]; }

    //! mixture density rho(c)
    Scalar density() const
    { return density_; }

    //! mixture (dynamic) viscosity nu(c)
    Scalar effectiveViscosity() const
    { return viscosity_; }

    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    const PrimaryVariables& priVars() const
    { return priVars_; }

private:
    PrimaryVariables priVars_;
    Scalar density_;
    Scalar viscosity_;
    Scalar extrusionFactor_;
};

} // end namespace Dumux

#endif
