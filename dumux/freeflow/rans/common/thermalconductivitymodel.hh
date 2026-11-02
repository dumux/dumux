// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::RANSThermalConductivityModel
 */
#ifndef DUMUX_RANS_COMMON_THERMAL_CONDUCTIVITY_MODEL_HH
#define DUMUX_RANS_COMMON_THERMAL_CONDUCTIVITY_MODEL_HH

#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Eddy thermal conductivity helper shared by the nonisothermal ("...NI") RANS mass
 *        TypeTags: \f$ \lambda_\text{t} = \mu_\text{t} c_\text{p} / \mathrm{Pr}_\text{t}
 *        \f$, related to the dynamic eddy viscosity \f$ \mu_\text{t} \f$ by the turbulent
 *        Prandtl number.
 *
 * \note `Properties::ThermalConductivityModel::effectiveThermalConductivity(volVars)` is called
 * from `Dumux::NavierStokesMassOnePVolumeVariables::update()` with `volVars` statically typed as
 * `NavierStokesMassOnePVolumeVariables<Traits>`, not as the derived RANS volume variables class,
 * because the former is the implementation type of `NavierStokesEnergyVolumeVariables`.
 * `dynamicEddyViscosity()`, which is only defined on the derived class, is therefore not
 * accessible from that call site.
 *
 * Hence, `Properties::ThermalConductivityModel` is set to the molecular-only
 * `effectiveThermalConductivity()` below, and the RANS mass volume variables call
 * `turbulentEffectiveThermalConductivity(volVars)` explicitly from their own `update()` with the
 * derived type, overwriting the protected `lambdaEff_` member of
 * `NavierStokesEnergyVolumeVariables` with the eddy-augmented value.
 */
struct RANSThermalConductivityModel
{
    //! Molecular-only effective thermal conductivity, to be used as
    //! Properties::ThermalConductivityModel (see class docs).
    template<class VolumeVariables>
    static auto effectiveThermalConductivity(const VolumeVariables& volVars)
    { return volVars.fluidThermalConductivity(); }

    //! The eddy-augmented effective thermal conductivity. Has to be called explicitly with the
    //! derived RANS volume variables type.
    template<class VolumeVariables>
    static auto turbulentEffectiveThermalConductivity(const VolumeVariables& volVars)
    {
        using FluidSystem = typename VolumeVariables::FluidSystem;
        using Scalar = std::decay_t<decltype(volVars.fluidThermalConductivity())>;
        static const Scalar turbulentPrandtlNumber = getParam<Scalar>("RANS.TurbulentPrandtlNumber", 1.0);

        return volVars.fluidThermalConductivity()
               + volVars.dynamicEddyViscosity()*FluidSystem::heatCapacity(volVars.fluidState(), 0)/turbulentPrandtlNumber;
    }
};

} // end namespace Dumux

#endif
