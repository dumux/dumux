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
