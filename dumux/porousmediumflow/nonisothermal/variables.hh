// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NIModel
 * \brief The energy part of the variables per local dof of porous-medium flow models
 */
#ifndef DUMUX_ENERGY_VARIABLES_HH
#define DUMUX_ENERGY_VARIABLES_HH

#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup NIModel
 * \brief The energy part of the variables per local dof of porous-medium flow models
 * \tparam Traits the traits of the volume variables of the model
 * \tparam Impl the variables class deriving from this class
 * \tparam enableEnergyBalance if the model solves an energy balance
 */
template<class Traits, class Impl, bool enableEnergyBalance = Traits::ModelTraits::enableEnergyBalance()>
class EnergyVariables;

/*!
 * \ingroup NIModel
 * \brief The energy part of the variables of isothermal models, with the temperature given by the spatial parameters
 */
template<class Traits, class Impl>
class EnergyVariables<Traits, Impl, false>
{
    using Scalar = typename Traits::PrimaryVariables::value_type;

public:
    using FluidState = typename Traits::FluidState;
    using SolidState = typename Traits::SolidState;
    using FluidSystem = typename Traits::FluidSystem;

    //! Sets the temperature of the spatial parameters in the fluid and the solid state
    template<class ElementSolution, class Problem, class ElementDiscretization, Concept::LocalDofIpData IpData>
    void updateTemperature(const ElementSolution& elemSol,
                           const Problem& problem,
                           const ElementDiscretization& elemDisc,
                           const IpData& ipData,
                           FluidState& fluidState,
                           SolidState& solidState)
    {
        const Scalar temperature = problem.spatialParams().temperature(elemDisc, ipData, elemSol);
        for (int phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
            fluidState.setTemperature(phaseIdx, temperature);
        solidState.setTemperature(temperature);
    }

    //! The solid energy parameters do not enter isothermal models
    template<class ElementSolution, class Problem, class ElementDiscretization, Concept::LocalDofIpData IpData>
    void updateSolidEnergyParams(const ElementSolution&, const Problem&, const ElementDiscretization&,
                                 const IpData&, SolidState&)
    {}

    //! The phase enthalpy is zero for isothermal models
    template<class ParameterCache>
    static Scalar enthalpy(const FluidState&, const ParameterCache&, const int)
    { return 0.0; }

    //! The effective thermal conductivity does not enter isothermal models
    void updateEffectiveThermalConductivity()
    {}
};

/*!
 * \ingroup NIModel
 * \brief The energy part of the variables of models with an energy balance in local thermal equilibrium,
 *        i.e. with one temperature for all phases and the solid
 *
 * The solid heat capacity, density and thermal conductivity are those of the solid system. For the
 * solid system SolidSystems::InertSolidPhase, the spatial parameters can overwrite them with
 * `solidHeatCapacity(elemDisc, ipData, elemSol, solidState)`, `solidDensity(...)` and
 * `solidThermalConductivity(...)` with the same arguments.
 */
template<class Traits, class Impl>
class EnergyVariables<Traits, Impl, true>
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using Indices = typename Traits::ModelTraits::Indices;
    using EffectiveThermalConductivityModel = typename Traits::EffectiveThermalConductivityModel;

    static_assert(Traits::ModelTraits::numEnergyEq() == 1,
                  "Local thermal non-equilibrium is not implemented for variables defined per local dof");

public:
    using FluidState = typename Traits::FluidState;
    using SolidState = typename Traits::SolidState;
    using FluidSystem = typename Traits::FluidSystem;
    using SolidSystem = typename Traits::SolidSystem;

    //! Sets the temperature of the primary variables in the fluid and the solid state
    template<class ElementSolution, class Problem, class ElementDiscretization, Concept::LocalDofIpData IpData>
    void updateTemperature(const ElementSolution& elemSol,
                           const Problem& problem,
                           const ElementDiscretization& elemDisc,
                           const IpData& ipData,
                           FluidState& fluidState,
                           SolidState& solidState)
    {
        const Scalar temperature = elemSol[ipData.localDofIndex()][Indices::temperatureIdx];
        for (int phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
            fluidState.setTemperature(phaseIdx, temperature);
        solidState.setTemperature(temperature);
    }

    //! Sets the heat capacity, density and thermal conductivity of the solid
    template<class ElementSolution, class Problem, class ElementDiscretization, Concept::LocalDofIpData IpData>
    void updateSolidEnergyParams(const ElementSolution& elemSol,
                                 const Problem& problem,
                                 const ElementDiscretization& elemDisc,
                                 const IpData& ipData,
                                 SolidState& solidState)
    {
        const auto& spatialParams = problem.spatialParams();

        if constexpr (requires { spatialParams.solidHeatCapacity(elemDisc, ipData, elemSol, solidState); })
        {
            static_assert(isInertSolidPhase_, "solidHeatCapacity can only be overwritten in the spatial params for the solid system InertSolidPhase");
            solidState.setHeatCapacity(spatialParams.solidHeatCapacity(elemDisc, ipData, elemSol, solidState));
        }
        else
            solidState.setHeatCapacity(SolidSystem::heatCapacity(solidState));

        if constexpr (requires { spatialParams.solidDensity(elemDisc, ipData, elemSol, solidState); })
        {
            static_assert(isInertSolidPhase_, "solidDensity can only be overwritten in the spatial params for the solid system InertSolidPhase");
            solidState.setDensity(spatialParams.solidDensity(elemDisc, ipData, elemSol, solidState));
        }
        else
            solidState.setDensity(SolidSystem::density(solidState));

        if constexpr (requires { spatialParams.solidThermalConductivity(elemDisc, ipData, elemSol, solidState); })
        {
            static_assert(isInertSolidPhase_, "solidThermalConductivity can only be overwritten in the spatial params for the solid system InertSolidPhase");
            solidState.setThermalConductivity(spatialParams.solidThermalConductivity(elemDisc, ipData, elemSol, solidState));
        }
        else
            solidState.setThermalConductivity(SolidSystem::thermalConductivity(solidState));
    }

    //! Updates the effective thermal conductivity of the fluid-solid aggregate
    void updateEffectiveThermalConductivity()
    { effectiveThermalConductivity_ = EffectiveThermalConductivityModel::effectiveThermalConductivity(asImp_()); }

    //! The specific enthalpy of a fluid phase
    template<class ParameterCache>
    static Scalar enthalpy(const FluidState& fluidState, const ParameterCache& paramCache, const int phaseIdx)
    { return FluidSystem::enthalpy(fluidState, paramCache, phaseIdx); }

    //! The specific internal energy of a fluid phase
    Scalar internalEnergy(const int phaseIdx) const
    { return asImp_().fluidState().internalEnergy(phaseIdx); }

    //! The specific enthalpy of a fluid phase
    Scalar enthalpy(const int phaseIdx) const
    { return asImp_().fluidState().enthalpy(phaseIdx); }

    //! The specific heat capacity of the solid
    Scalar solidHeatCapacity() const
    { return asImp_().solidState().heatCapacity(); }

    //! The mass density of the solid
    Scalar solidDensity() const
    { return asImp_().solidState().density(); }

    //! The thermal conductivity of the solid
    Scalar solidThermalConductivity() const
    { return asImp_().solidState().thermalConductivity(); }

    //! The thermal conductivity of a fluid phase
    Scalar fluidThermalConductivity(const int phaseIdx) const
    { return FluidSystem::thermalConductivity(asImp_().fluidState(), phaseIdx); }

    //! The effective thermal conductivity of the fluid-solid aggregate
    Scalar effectiveThermalConductivity() const
    { return effectiveThermalConductivity_; }

private:
    static constexpr bool isInertSolidPhase_ = Dumux::Detail::isInertSolidPhase<SolidSystem>::value;

    const Impl& asImp_() const { return *static_cast<const Impl*>(this); }

    Scalar effectiveThermalConductivity_;
};

} // end namespace Dumux::Experimental

#endif
