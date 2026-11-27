// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \brief Helper struct defining the advective fluxes of the one-equation (Spalart-Allmaras)
 *        turbulence closure on top of the single-phase Navier-Stokes mass model, for use
 *        with the boundary flux helpers in dumux/freeflow/navierstokes/scalarfluxhelper.hh.
 */
#ifndef DUMUX_RANS_ONEEQ_MASS_ADVECTIVE_FLUX_HH
#define DUMUX_RANS_ONEEQ_MASS_ADVECTIVE_FLUX_HH

// This include provides the primary template (including its default template argument), which
// must not be redeclared here since default template arguments cannot be repeated.
#include <dumux/freeflow/navierstokes/mass/1p/advectiveflux.hh>
#include <dumux/freeflow/navierstokes/energy/model.hh>

namespace Dumux {

#ifndef DOXYGEN
struct OneEqMassModelTraits;
#endif

/*!
 * \ingroup FreeflowModels
 * \brief Helper struct defining the advective fluxes of the one-equation turbulence closure
 *        on top of the single-phase Navier-Stokes mass model: the mass flux of the
 *        pressure equation (density upwind) plus the advective flux of the ν̃ transport
 *        equation (ρν̃ upwind).
 */
template<class T>
struct AdvectiveFlux<OneEqMassModelTraits, T>
{
    template<class NumEqVector, class UpwindFunction>
    static void addAdvectiveFlux(NumEqVector& flux,
                                 const UpwindFunction& upwind)
    {
        using ModelTraits = T;

        flux[ModelTraits::Indices::conti0EqIdx] += upwind([](const auto& volVars) { return volVars.density(); });
        flux[ModelTraits::Indices::viscosityTildeEqIdx]
            += upwind([](const auto& volVars) { return volVars.density()*volVars.viscosityTilde(); });
    }
};

//! Use the same advective flux for the non-isothermal model (heat fluxes are added separately).
//! Note: OneEqMassLocalResidual overwrites the viscosityTildeEqIdx entry of the flux, so only
//! the conti0EqIdx entry is relevant here.
template<>
struct AdvectiveFlux<NavierStokesEnergyModelTraits<OneEqMassModelTraits>>
: public AdvectiveFlux<OneEqMassModelTraits>
{};

} // end namespace Dumux

#endif
