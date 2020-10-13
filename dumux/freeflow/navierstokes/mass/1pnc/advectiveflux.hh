// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
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
        static constexpr bool useTotalMoleOrMassBalance = replaceCompEqIdx < numComponents;

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
        if constexpr(useTotalMoleOrMassBalance)
        {
            auto upwindTerm = [&]()
            {
                if constexpr (useMoles)
                    return [](const auto& volVars) { return volVars.molarDensity(); };
                else
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
