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
 * \ingroup FreeflowNIModel
 * \copydoc Dumux::FreeFlowEnergyLocalResidual
 */
#ifndef DUMUX_NAVIERSTOKES_ENERGY_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_ENERGY_LOCAL_RESIDUAL_HH

#include <dumux/discretization/method.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {


/*!
 * \ingroup FreeflowNIModel
 * \brief DOCME
 */
template<class ModelTraits>
class NavierStokesEnergyLocalResidual
{
    using Indices = typename ModelTraits::Indices;
    static constexpr bool enableEnergyBalance = ModelTraits::enableEnergyBalance();
public:

    //! The energy storage in the fluid phase
    template <class NumEqVector, class VolumeVariables>
    static void fluidPhaseStorage(NumEqVector& storage,
                                  const VolumeVariables& volVars)
    {
        if constexpr (enableEnergyBalance)
            storage[Indices::energyEqIdx] = volVars.density() * volVars.internalEnergy();
    }

    //! brief The advective phase energy flux
    template <class NumEqVector, class FluxVariables>
    static void heatConvectionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars)
    {
        if constexpr (enableEnergyBalance)
        {
            auto upwindTerm = [](const auto& volVars)
            { return volVars.density()*volVars.enthalpy(); };

            flux[Indices::energyEqIdx] += fluxVars.advectiveFlux(upwindTerm);
        }
    }

    //! brief The conductive phase energy flux
    template <class NumEqVector, class FluxVariables>
    static void heatConductionFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars)
    {
        if constexpr (enableEnergyBalance)
            flux[Indices::energyEqIdx] += fluxVars.heatConductionFlux();
    }

    // TODO diffusive enthalpy flux
};


} // end namespace Dumux

#endif
