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
 * \brief Eddy thermal conductivity helper shared by every nonisothermal ("...NI") RANS mass
 *        TypeTag: \f$ \lambda_\text{t} = \nu_\text{t} \varrho c_\text{p} / \mathrm{Pr}_\text{t}
 *        \f$, related to the eddy viscosity by the turbulent Prandtl number - the direct
 *        analogue of the deleted releases/3.10:dumux/freeflow/rans/volumevariables.hh's
 *        `calculateEddyThermalConductivity()`.
 *
 * \note `Properties::ThermalConductivityModel::effectiveThermalConductivity(volVars)` (called
 * from `Dumux::NavierStokesEnergyVolumeVariables::updateEffectiveThermalConductivity()`, in turn
 * called from *inside* `Dumux::NavierStokesMassOnePVolumeVariables::update()`) is invoked with
 * `volVars` statically typed as `NavierStokesMassOnePVolumeVariables<Traits>` - **not** the
 * further-derived RANS volVars class (e.g. `Dumux::KOmegaMassVolumeVariables<Traits>`) - because
 * `NavierStokesMassOnePVolumeVariables` hardcodes itself (not a passed-through template
 * parameter) as the CRTP `Impl` type of `NavierStokesEnergyVolumeVariables<Traits, Impl>`. So
 * `dynamicEddyViscosity()` (defined only on the derived RANS class) is not reachable from that
 * generic call site - unlike the CRTP hooks used elsewhere in this project (e.g.
 * `RANSMomentumProblem::isOnWall()`), there is no existing "goes-through-`asImp_()`-correctly"
 * forwarding method to reuse here, since this one lives in shared DuMux framework code, not in
 * project code.
 *
 * The fix used throughout this project's RANS mass volume-variables classes: keep
 * `Properties::ThermalConductivityModel` set to a harmless, always-compiling passthrough
 * (`effectiveThermalConductivity()` below, molecular-only, safe against the base-sliced type),
 * and instead call `turbulentEffectiveThermalConductivity(volVars)` **explicitly**, with the
 * correctly, fully-derived `volVars` type, from each RANS model's own
 * `update()`/`updateRANSProperties()` - overwriting the protected `lambdaEff_` member (inherited,
 * and thus still accessible, from `NavierStokesEnergyVolumeVariables`) with the eddy-augmented
 * value. See e.g. `dumux/freeflow/rans/komega/volumevariables.hh`.
 */
struct RANSThermalConductivityModel
{
    //! Safe to wire as Properties::ThermalConductivityModel: molecular-only, never touches
    //! dynamicEddyViscosity() (see class docs for why it structurally cannot).
    template<class VolumeVariables>
    static auto effectiveThermalConductivity(const VolumeVariables& volVars)
    { return volVars.fluidThermalConductivity(); }

    //! The actual eddy-augmented effective thermal conductivity - call this explicitly, with the
    //! fully-derived RANS volVars type, from that class's own update()/updateRANSProperties().
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
