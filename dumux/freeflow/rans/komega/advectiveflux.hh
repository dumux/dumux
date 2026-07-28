// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \brief Helper struct defining the advective fluxes of the two-equation k-omega turbulence
 *        closure fused onto the single-phase Navier-Stokes mass model, for use with
 *        dumux/freeflow/navierstokes/scalarfluxhelper.hh's boundary flux helpers (the outflow
 *        treatment at the test's outlet, see test/freeflow/rans/komega/channel).
 */
#ifndef DUMUX_RANS_KOMEGA_MASS_ADVECTIVE_FLUX_HH
#define DUMUX_RANS_KOMEGA_MASS_ADVECTIVE_FLUX_HH

// Pulls in the primary template (with its default second argument) exactly once - see
// dumux/freeflow/rans/oneeq/advectiveflux.hh for why this must not be redeclared here.
#include <dumux/freeflow/navierstokes/mass/1p/advectiveflux.hh>
#include <dumux/freeflow/navierstokes/energy/model.hh>

namespace Dumux {

#ifndef DOXYGEN
struct KOmegaMassModelTraits;
#endif

/*!
 * \ingroup FreeflowModels
 * \brief Helper struct defining the advective fluxes of the k-omega turbulence closure fused
 *        onto the single-phase Navier-Stokes mass model: the pressure equation's mass flux
 *        (density upwind, as in the plain 1p model) plus k's and omega's advective fluxes
 *        (rho*k, rho*omega upwind).
 */
template<class T>
struct AdvectiveFlux<KOmegaMassModelTraits, T>
{
    template<class NumEqVector, class UpwindFunction>
    static void addAdvectiveFlux(NumEqVector& flux,
                                 const UpwindFunction& upwind)
    {
        using ModelTraits = T;

        flux[ModelTraits::Indices::conti0EqIdx] += upwind([](const auto& volVars) { return volVars.density(); });
        flux[ModelTraits::Indices::turbulentKineticEnergyEqIdx]
            += upwind([](const auto& volVars) { return volVars.density()*volVars.turbulentKineticEnergy(); });
        flux[ModelTraits::Indices::dissipationEqIdx]
            += upwind([](const auto& volVars) { return volVars.density()*volVars.dissipation(); });
    }
};

//! Use the same advective flux for the non-isothermal model (heat fluxes are added separately
//! by NavierStokesMassOnePFluxVariables::flux(), see dumux/freeflow/rans/komega/model.hh's NI
//! TypeTag) - mirrors dumux/freeflow/navierstokes/mass/1p/advectiveflux.hh's own NI passthrough.
template<>
struct AdvectiveFlux<NavierStokesEnergyModelTraits<KOmegaMassModelTraits>>
: public AdvectiveFlux<KOmegaMassModelTraits>
{};

} // end namespace Dumux

#endif
