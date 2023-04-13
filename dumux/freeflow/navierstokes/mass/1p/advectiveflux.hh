// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
 /*!
  * \file
  * \ingroup NavierStokesModel
  * \brief Helper struct defining the advective fluxes of the single-phase flow
  *        Navier-Stokes mass model
  */
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_MASS_1P_ADVECTIVE_FLUX_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_MASS_1P_ADVECTIVE_FLUX_HH

namespace Dumux {

#ifndef DOXYGEN
// forward declare
struct NavierStokesMassOnePModelTraits;

template<class IsothermalTraits>
struct NavierStokesEnergyModelTraits;

template<class ModelTraits, class T = ModelTraits>
struct AdvectiveFlux;
#endif

/*!
 * \ingroup NavierStokesModel
 * \brief Helper struct defining the advective fluxes of the single-phase flow
 *        Navier-Stokes mass model
 */
template<class T>
struct AdvectiveFlux<NavierStokesMassOnePModelTraits, T>
{
    template<class NumEqVector, class UpwindFunction>
    static void addAdvectiveFlux(NumEqVector& flux,
                                 const UpwindFunction& upwind)
    {
        using ModelTraits = T;

        // get equation index
        const auto eqIdx = ModelTraits::Indices::conti0EqIdx;
        flux[eqIdx] += upwind([](const auto& volVars) { return volVars.density(); });
    }
};

// use the same mass flux for the non-isothermal model (heat fluxes are added separately)
template<>
struct AdvectiveFlux<NavierStokesEnergyModelTraits<NavierStokesMassOnePModelTraits>>
: public AdvectiveFlux<NavierStokesMassOnePModelTraits>
{};

} // end namespace Dumux

#endif
