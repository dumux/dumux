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
 *        turbulence closure fused onto the single-phase Navier-Stokes mass model, for use
 *        with dumux/freeflow/navierstokes/scalarfluxhelper.hh's boundary flux helpers
 *        (e.g. the outflow treatment at the test's outlet, see test/freeflow/rans/oneeq/channel).
 */
#ifndef DUMUX_RANS_ONEEQ_MASS_ADVECTIVE_FLUX_HH
#define DUMUX_RANS_ONEEQ_MASS_ADVECTIVE_FLUX_HH

// Pulls in the primary template (with its default second argument) exactly once, so this
// header does not need to (and must not, see below) redeclare it itself: repeating a default
// template argument for an already-visible declaration is a hard error in C++.
#include <dumux/freeflow/navierstokes/mass/1p/advectiveflux.hh>

namespace Dumux {

#ifndef DOXYGEN
struct OneEqMassModelTraits;
#endif

/*!
 * \ingroup FreeflowModels
 * \brief Helper struct defining the advective fluxes of the one-equation turbulence closure
 *        fused onto the single-phase Navier-Stokes mass model: the pressure equation's mass
 *        flux (density upwind, as in the plain 1p model) plus the ν̃ transport equation's
 *        advective flux (ρν̃ upwind).
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

} // end namespace Dumux

#endif
