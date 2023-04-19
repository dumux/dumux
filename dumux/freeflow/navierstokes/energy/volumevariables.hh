// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowNIModel
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities.
 */

#ifndef DUMUX_FREEFLOW_NAVIER_STOKES_ENERGY_VOLUME_VARIABLES_HH
#define DUMUX_FREEFLOW_NAVIER_STOKES_ENERGY_VOLUME_VARIABLES_HH

#include <type_traits>
#include <dune/common/std/type_traits.hh>


namespace Dumux {

namespace Detail {

struct EmptyFreeFlowHeatCondType {};

template<bool enableEnergyBalance, class Traits>
struct FreeFlowHeatCondType
{
    using type = EmptyFreeFlowHeatCondType;
};

template<class Traits>
struct FreeFlowHeatCondType<true, Traits>
{
    using type = typename Traits::HeatConductionType;
};

} // end namespace Detail

/*!
 * \ingroup FreeflowNIModel
 * \brief The isothermal base class
 */
template<class Traits, class Impl>
class NavierStokesEnergyVolumeVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    static constexpr bool enableEnergyBalance = Traits::ModelTraits::enableEnergyBalance();

public:
    using FluidState = typename Traits::FluidState;
    using FluidSystem = typename Traits::FluidSystem;
    using HeatConductionType = typename Detail::FreeFlowHeatCondType<enableEnergyBalance, Traits>::type;

    /*!
    * \brief Returns the temperature at a given sub-control volume
    *
    * \param elemSol A vector containing all primary variables connected to the element
    * \param problem The object specifying the problem which ought to
    *                be simulated
    * \param element An element which contains part of the control volume
    * \param scv The sub-control volume
    */
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    Scalar getTemperature(const ElementSolution& elemSol,
                          const Problem& problem,
                          const Element& element,
                          const SubControlVolume& scv) const
    {
        if constexpr (enableEnergyBalance)
            return elemSol[scv.localDofIndex()][Traits::ModelTraits::Indices::temperatureIdx];
        else
            return problem.spatialParams().temperature(element, scv, elemSol);
    }


    //! The effective thermal conductivity is zero for isothermal models
    void updateEffectiveThermalConductivity()
    {
        if constexpr (enableEnergyBalance)
            lambdaEff_ = Traits::EffectiveThermalConductivityModel::effectiveThermalConductivity(asImp_());
    }

    /*!
     * \brief Returns the total internal energy of a phase in the
     *        sub-control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar internalEnergy(const int phaseIdx = 0) const
    {
        if constexpr (enableEnergyBalance)
            return asImp_().fluidState().internalEnergy(0);
        else
            return 0.0;
    }

    /*!
     * \brief Returns the total enthalpy of a phase in the sub-control
     *        volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar enthalpy(const int phaseIdx = 0) const
    {
        if constexpr (enableEnergyBalance)
            return asImp_().fluidState().enthalpy(0);
        else
            return 0.0;
    }

    /*!
     * \brief Returns the temperature of a fluid phase assuming thermal nonequilibrium
     *        the sub-control volume.
     * \param phaseIdx The local index of the phases
     */
    Scalar temperatureFluid(const int phaseIdx = 0) const
    { return asImp_().fluidState().temperature(0); }


    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m*K)]}\f$
     *        of a fluid phase in the sub-control volume.
     */
    Scalar fluidThermalConductivity(const int phaseIdx = 0) const
    { return FluidSystem::thermalConductivity(asImp_().fluidState(), 0); }

    /*!
     * \brief Returns the effective thermal conductivity \f$\mathrm{[W/(m*K)]}\f$ in
     *        the sub-control volume. Specific to equilibirum models (case fullThermalEquilibrium).
     */
    Scalar effectiveThermalConductivity(const int phaseIdx = 0) const
    { return lambdaEff_; }


    //! The phase enthalpy is zero for isothermal models
    //! This is needed for completing the fluid state
    template<class ParameterCache>
    static Scalar enthalpy(const FluidState& fluidState,
                           const ParameterCache& paramCache)
    {
        if constexpr (enableEnergyBalance)
            return FluidSystem::enthalpy(fluidState, paramCache, 0);
        else
            return 0.0;
    }

protected:
    Scalar lambdaEff_;
    const Impl &asImp_() const { return *static_cast<const Impl*>(this); }
    Impl &asImp_() { return *static_cast<Impl*>(this); }
};

} // end namespace Dumux

#endif
