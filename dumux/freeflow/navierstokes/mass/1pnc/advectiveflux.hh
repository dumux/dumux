// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
 /*!
  * \file
  * \ingroup NavierStokesModel
  * \brief Helper struct defining the advective fluxes of the single-phase flow multicomponent
  *        Navier-Stokes mass model
  */
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_MASS_1PNC_ADVECTIVE_FLUX_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_MASS_1PNC_ADVECTIVE_FLUX_HH

namespace Dumux {

#ifndef DOXYGEN
// forward declare
template<int nComp, bool useM, int repCompEqIdx>
struct NavierStokesMassOnePNCModelTraits;

template<class IsothermalTraits>
struct NavierStokesEnergyModelTraits;

template<class ModelTraits>
struct AdvectiveFlux;
#endif


/*!
 * \ingroup NavierStokesModel
 * \brief Helper struct defining the advective fluxes of the single-phase flow multicomponent
 *        Navier-Stokes mass model
 */
template<int nComp, bool useM, int repCompEqIdx>
struct AdvectiveFlux<NavierStokesMassOnePNCModelTraits<nComp, useM, repCompEqIdx>>
{
    template<class NumEqVector, class UpwindFunction>
    static void addAdvectiveFlux(NumEqVector& flux,
                                 const UpwindFunction& upwind)
    {
        using ModelTraits = NavierStokesMassOnePNCModelTraits<nComp, useM, repCompEqIdx>;
        static constexpr bool useMoles = ModelTraits::useMoles();
        static constexpr auto numComponents = ModelTraits::numFluidComponents();
        static constexpr auto replaceCompEqIdx = ModelTraits::replaceCompEqIdx();
        static constexpr bool useTotalMassBalance = replaceCompEqIdx < numComponents;

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            // get equation index
            const auto eqIdx = ModelTraits::Indices::conti0EqIdx + compIdx;

            if (eqIdx == replaceCompEqIdx)
                continue;

            auto upwindTerm = [&]()
            {
                if constexpr (useMoles)
                    return [compIdx](const auto& volVars) { return volVars.molarDensity()*volVars.moleFraction(compIdx); };
                else
                    return [compIdx](const auto& volVars) { return volVars.density()*volVars.massFraction(compIdx); };
            }();

            flux[eqIdx] += upwind(upwindTerm);
        }

        // in case one balance is substituted by the total mole balance
        if constexpr(useTotalMassBalance)
        {
            auto upwindTerm = [&]()
            {
                return [](const auto& volVars) { return volVars.density(); };
            }();

            flux[replaceCompEqIdx] += upwind(upwindTerm);
        }
    }
};

// use the same mass flux for the non-isothermal model (heat fluxes are added separately)
template<int nComp, bool useM, int repCompEqIdx>
struct AdvectiveFlux<NavierStokesEnergyModelTraits<NavierStokesMassOnePNCModelTraits<nComp, useM, repCompEqIdx>>>
: public AdvectiveFlux<NavierStokesMassOnePNCModelTraits<nComp, useM, repCompEqIdx>>
{};

} // end namespace Dumux

#endif
