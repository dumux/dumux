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
 * \ingroup Fluidsystems
 * \brief @copybrief Dumux::FluidSystems::OnePAdapter
 */
#ifndef DUMUX_FLUIDSYTEMS_ONEP_ADAPTER_HH
#define DUMUX_FLUIDSYTEMS_ONEP_ADAPTER_HH

#include <cassert>

#include <dune/common/exceptions.hh>

#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/fluidstates/adapter.hh>

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 * \brief An adapter for multi-phase fluid systems to be used with (compositional) one-phase models
 * \tparam MPFluidSystem the multi-phase fluid system to be adapted
 * \tparam phase the index of the phase we choose from the multi-phase fluid system
 */
template <class MPFluidSystem, int phase = 0>
class OnePAdapter
: public Base<typename MPFluidSystem::Scalar, OnePAdapter<MPFluidSystem, phase>>
{
    using ThisType = OnePAdapter<MPFluidSystem, phase>;
    using Base = Dumux::FluidSystems::Base<typename MPFluidSystem::Scalar, ThisType>;

    static_assert(phase < MPFluidSystem::numPhases, "Phase does not exist in multi-phase fluidsystem!");

    struct AdapterPolicy
    {
        using FluidSystem = MPFluidSystem;

        // the phase index is always zero, other phases than the chosen phase should never be called
        static int phaseIdx(int mpFluidPhaseIdx)
        { return 0; }

        // the main component is currently excepted to have the same index as it's phase
        // (see Fluidsystems::Base::getMainComponent for more information)
        // so we swap the main component with the first component
        // this mapping works in both ways since we are only swapping components
        static constexpr int compIdx(int compIdx)
        {
            if (compIdx == 0)
                return phase;
            else if (compIdx == phase)
                return 0;
            else
                return compIdx;
        }
    };

    template<class FluidState>
    static auto adaptFluidState(const FluidState& fluidState)
    { return FluidStateAdapter<FluidState, AdapterPolicy>(fluidState); }

public:
    using Scalar = typename Base::Scalar;
    using ParameterCache = NullParameterCache;

    //! export the wrapped MultiPhaseFluidSystem type
    using MultiPhaseFluidSystem = MPFluidSystem;
    //! the index of the phase we choose from the multi-phase fluid system
    static constexpr int multiphaseFluidsystemPhaseIdx = phase;

    //! number of phases in the fluid system
    static constexpr int numPhases = 1;
    //! for compositional models, the number of components has to be the same as in the multi-phase fluid system as the composition needs to be defined,
    //! while for non-compositional models, the number of components must equal the number of phases (1 in this case)
    static constexpr int numComponents = MultiPhaseFluidSystem::isMiscible() ? MultiPhaseFluidSystem::numComponents : numPhases;
    //! number of components has to be the same as in the multi-phase fluid system as the composition needs to be defined
    static constexpr int phase0Idx = 0; //!< index of the only phase

    //! convert a component index of the multi-phase component index to the actual component index
    static constexpr int compIdx(int multiPhaseFluidSystemCompIdx)
    { return AdapterPolicy::compIdx(multiPhaseFluidSystemCompIdx); }

    /*!
     * \brief Initialize the fluid system's static parameters generically
     */
    template<class ...Args>
    static void init(Args&&... args)
    { MultiPhaseFluidSystem::init(std::forward<Args>(args)...); }

    /****************************************
     * Fluid phase related static parameters
     ****************************************/
    /*!
     * \brief Return the human readable name of a fluid phase
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx = 0)
    { return MultiPhaseFluidSystem::phaseName(phase); }

    /*!
     * \brief A human readable name for the component.
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    { return MultiPhaseFluidSystem::componentName(AdapterPolicy::compIdx(compIdx)); }

    /*!
     * \brief A human readable name for the component.
     */
    static std::string name()
    { return MultiPhaseFluidSystem::phaseName(phase); }

    /*!
     * \brief There is only one phase, so not mass transfer between phases can occur
     */
    static constexpr bool isMiscible()
    { return false; }

    /*!
     * \brief Returns whether the fluid is gaseous
     */
    static constexpr bool isGas(int phaseIdx = 0)
    { return MultiPhaseFluidSystem::isGas(phase); }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are independent on the fluid composition. This assumption is true
     * if only a single component is involved. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isIdealMixture(int phaseIdx = 0)
    { return MultiPhaseFluidSystem::isIdealMixture(phase); }

    /*!
     * \brief Returns true if the fluid is assumed to be compressible
     */
    static constexpr bool isCompressible(int phaseIdx = 0)
    { return MultiPhaseFluidSystem::isCompressible(phase); }

    /*!
     * \brief Returns true if the fluid viscosity is constant
     */
    static constexpr bool viscosityIsConstant(int phaseIdx = 0)
    { return MultiPhaseFluidSystem::viscosityIsConstant(phase); }

    /*!
     * \brief Returns true if the fluid is assumed to be an ideal gas
     */
    static constexpr bool isIdealGas(int phaseIdx = 0)
    { return MultiPhaseFluidSystem::isIdealGas(phase); }

    /*!
     * \brief The mass in \f$\mathrm{[kg]}\f$ of one mole of the component.
     */
    static Scalar molarMass(int compIdx)
    {  return MultiPhaseFluidSystem::molarMass(AdapterPolicy::compIdx(compIdx)); }

    using Base::density;
    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure and temperature.
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState, int phaseIdx = 0)
    {
        assert(phaseIdx == 0);
        return MultiPhaseFluidSystem::density(adaptFluidState(fluidState), phase);
    }

    using Base::molarDensity;
    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density is defined by the
     * mass density \f$\rho_\alpha\f$ and the main component molar mass \f$M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{M_\alpha} \;.\f]
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, int phaseIdx = 0)
    {
        assert(phaseIdx == 0);
        return MultiPhaseFluidSystem::molarDensity(adaptFluidState(fluidState), phase);
    }

    using Base::enthalpy;
    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ the pure component as a liquid.
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState, int phaseIdx = 0)
    {
        assert(phaseIdx == 0);
        return MultiPhaseFluidSystem::enthalpy(adaptFluidState(fluidState), phase);
    }

    /*!
     * \brief Returns the specific enthalpy \f$\mathrm{[J/kg]}\f$ of a component in a specific phase
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     *
     */
    template <class FluidState>
    static Scalar componentEnthalpy(const FluidState &fluidState,
                                    int phaseIdx,
                                    int compIdx)
    {
        assert(phaseIdx == 0);
        return MultiPhaseFluidSystem::componentEnthalpy(adaptFluidState(fluidState), phase,
                                                        AdapterPolicy::compIdx(compIdx));
    }

    using Base::viscosity;
    /*!
     * \brief The dynamic liquid viscosity \f$\mathrm{[N/m^3*s]}\f$ of the pure component.
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState, int phaseIdx = 0)
    {
        assert(phaseIdx == 0);
        return MultiPhaseFluidSystem::viscosity(adaptFluidState(fluidState), phase);
    }

    using Base::fugacityCoefficient;
    /*!
     * \copybrief Base::fugacityCoefficient
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(phaseIdx == 0);
        return MultiPhaseFluidSystem::fugacityCoefficient(adaptFluidState(fluidState), phase,
                                                          AdapterPolicy::compIdx(compIdx));
    }

    using Base::diffusionCoefficient;
    /*!
     * \copybrief Base::diffusionCoefficient
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        assert(phaseIdx == 0);
        return MultiPhaseFluidSystem::diffusionCoefficient(adaptFluidState(fluidState), phase,
                                                           AdapterPolicy::compIdx(compIdx));
    }

    using Base::binaryDiffusionCoefficient;
    /*!
     * \copybrief Base::binaryDiffusionCoefficient
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIIdx The index of the component to consider
     * \param compJIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        assert(phaseIdx == 0);
        return MultiPhaseFluidSystem::binaryDiffusionCoefficient(adaptFluidState(fluidState), phase,
                                                                 AdapterPolicy::compIdx(compIIdx),
                                                                 AdapterPolicy::compIdx(compJIdx));
    }

    using Base::thermalConductivity;
    /*!
     * \brief Thermal conductivity of the fluid \f$\mathrm{[W/(m K)]}\f$.
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      int phaseIdx = 0)
    {
        assert(phaseIdx == 0);
        return MultiPhaseFluidSystem::thermalConductivity(adaptFluidState(fluidState), phase);
    }

    using Base::heatCapacity;
    /*!
     * \brief Specific isobaric heat capacity of the fluid \f$\mathrm{[J/(kg K)]}\f$.
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx = 0)
    {
        assert(phaseIdx == 0);
        return MultiPhaseFluidSystem::heatCapacity(adaptFluidState(fluidState), phase);
    }
};

} // namespace FluidSystems
} // namespace

#endif
