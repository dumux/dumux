// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FluidSystems
 * \brief @copybrief Dumux::FluidSystems::GasPhaseTwoC
 */
#ifndef DUMUX_GAS_PHASE_HH
#define DUMUX_GAS_PHASE_HH

#include <cassert>
#include <limits>

#include <dune/common/exceptions.hh>

#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/components/componenttraits.hh>
#include <dumux/io/name.hh>

#include <dumux/material/idealgas.hh>

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup FluidSystems
 * \brief Policy for the GasPhaseTwoC fluid system
 */
template<bool fastButSimplifiedRelations = false>
struct GasPhaseTwoCDefaultPolicy
{
    static constexpr  bool useIdealGasDensity() { return fastButSimplifiedRelations; }
};

/*!
 * \ingroup FluidSystems
 * \brief A gaseous phase consisting of two components
 */
template <class Scalar, class MainComponent, class SecondComponent, class Policy = GasPhaseTwoCDefaultPolicy<> >
class GasPhaseTwoC
: public Base<Scalar, GasPhaseTwoC<Scalar, MainComponent, SecondComponent, Policy> >
{
    using ThisType = GasPhaseTwoC<Scalar, MainComponent, SecondComponent, Policy>;

    static_assert(ComponentTraits<MainComponent>::hasGasState, "The component does not implement a gas state!");
    static_assert(ComponentTraits<SecondComponent>::hasGasState, "The component does not implement a gas state!");

    // convenience aliases using declarations
    using IdealGas = Dumux::IdealGas<Scalar>;

public:
    using ParameterCache = NullParameterCache;

    static constexpr int numPhases = 1;  //!< Number of phases in the fluid system
    static constexpr int numComponents = 2; //!< Number of components in the fluid system

    static constexpr int phase0Idx = 0; //!< index of the only phase
    static constexpr int comp0Idx = 0; //!< index of the first component
    static constexpr int comp1Idx = 1; //!< index of the second component

    /*!
     * \brief Initialize the fluid system's static parameters generically
     */
    static void init()
    { }

    /****************************************
     * Fluid phase related static parameters
     ****************************************/
    /*!
     * \brief Return the human readable name of a fluid phase
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx = 0)
    { return IOName::gaseousPhase(); }

    /*!
     * \brief A human readable name for the component.
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        switch (compIdx)
        {
            case comp0Idx: return MainComponent::name();
            case comp1Idx: return SecondComponent::name();
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief A human readable name for the component.
     */
    static std::string name()
    { return "GasPhaseTwoC"; }

    /*!
     * \brief There is only one phase, so not mass transfer between phases can occur
     */
    static constexpr bool isMiscible()
    { return false; }

    /*!
     * \brief Returns whether the fluid is gaseous
     */
    static constexpr bool isGas(int phaseIdx = 0)
    { return true; }

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
    { return true; }

    /*!
     * \brief Returns true if the fluid is assumed to be compressible
     */
    static constexpr bool isCompressible(int phaseIdx)
    {
        return true;
    }

    /*!
     * \brief Returns true if the fluid is assumed to be an ideal gas
     */
    static constexpr bool isIdealGas(int phaseIdx,int compIdx)
    {
        switch(compIdx)
        {
            case comp0Idx: return MainComponent::gasIsIdeal();
            case comp1Idx: return SecondComponent::gasIsIdeal();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the molar mass of a component in \f$\mathrm{[kg/mol]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar molarMass(int compIdx)
    {
        static const Scalar M[] = {
            MainComponent::molarMass(),
            SecondComponent::molarMass(),
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return M[compIdx];
    }

    /*!
     * \brief Critical temperature of a component \f$\mathrm{[K]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalTemperature(int compIdx )
    {
        switch (compIdx)
        {
            case comp0Idx: return MainComponent::criticalTemperature();
            case comp1Idx: return SecondComponent::criticalTemperature();
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Critical pressure of a component \f$\mathrm{[Pa]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalPressure(int compIdx)
    {
        switch (compIdx)
        {
            case comp0Idx: return MainComponent::criticalPressure();
            case comp1Idx: return SecondComponent::criticalPressure();
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Returns the temperature in \f$\mathrm{[K]}\f$ at the component's triple point.
     */
    static Scalar tripleTemperature(int compIdx)
    {
        switch (compIdx)
        {
            case comp0Idx: return MainComponent::tripleTemperature();
            case comp1Idx: return SecondComponent::tripleTemperature();
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Returns the pressure in \f$\mathrm{[Pa]}\f$ at the component's triple point.
     */
    static Scalar triplePressure(int compIdx)
    {
        switch (compIdx)
        {
            case comp0Idx: return MainComponent::triplePressure();
            case comp1Idx: return SecondComponent::triplePressure();
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief The vapor pressure in \f$\mathrm{[Pa]}\f$ of the component at a given
     *        temperature.
     * \param T temperature \f$\mathrm{[K]}\f$
     */
    static Scalar vaporPressure(Scalar T, int compIdx)
    {
        switch (compIdx)
        {
            case comp0Idx: return MainComponent::vaporPressure();
            case comp1Idx: return SecondComponent::vaporPressure();
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure and temperature.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar density(Scalar temperature, Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "Density");
    }

    using Base<Scalar, ThisType>::density;
    //! \copydoc Base<Scalar,ThisType>::density(const FluidState&,int)
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        // for the gas phase assume an ideal gas
        using std::max;
        if (Policy::useIdealGasDensity())
            return IdealGas::molarDensity(T, p) * fluidState.averageMolarMass(0);

        // assume ideal mixture: steam, nitrogen and oxygen don't "see" each other
        else
            return MainComponent::gasDensity(T, fluidState.partialPressure(phaseIdx, comp0Idx))
            + SecondComponent::gasDensity(T, fluidState.partialPressure(phaseIdx, comp1Idx));
    }

    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density is defined by the
     * mass density \f$\rho_\alpha\f$ and the component molar mass \f$M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{M_\alpha} \;.\f]
     * \param temperature The temperature at which to evaluate the molar density
     * \param pressure The pressure at which to evaluate the molar density
     */
    static Scalar molarDensity(Scalar temperature, Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "Molar Density");
    }

    using Base<Scalar, ThisType>::molarDensity;
    //! \copydoc Base<Scalar,ThisType>::molarDensity(const FluidState&,int)
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState,
                               const int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        if (Policy::useIdealGasDensity())
        {   //assume ideal gas
            return IdealGas::molarDensity(T,p);
        }

        return MainComponent::gasMolarDensity(T, fluidState.partialPressure(phaseIdx, comp0Idx))
             + SecondComponent::gasMolarDensity(T, fluidState.partialPressure(phaseIdx, comp1Idx));
    }

    /*!
     * \brief The pressure \f$\mathrm{[Pa]}\f$ of the component at a given density and temperature.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param density The given density \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar pressure(Scalar temperature, Scalar density)
    {
        DUNE_THROW(Dune::NotImplemented, "Pressure");
    }

    /*!
     * \brief Specific enthalpy \f$\mathrm{[J/kg]}\f$ of the pure component as a gas.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static const Scalar enthalpy(Scalar temperature, Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "Enthalpy");
    }

    using Base<Scalar, ThisType>::enthalpy;
    //! \copydoc Base<Scalar,ThisType>::enthalpy(const FluidState&,int)
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           const int phaseIdx)
    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        return MainComponent::gasEnthalpy(T,p)*fluidState.massFraction(phaseIdx,comp0Idx)
            + SecondComponent::gasEnthalpy(T,p)*fluidState.massFraction(phaseIdx,comp1Idx);
    }

    /*!
     * \brief Returns the specific enthalpy \f$\mathrm{[J/kg]}\f$ of a component in a specific phase
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param componentIdx The index of the component to consider
     *
     */
    template <class FluidState>
    static Scalar componentEnthalpy(const FluidState &fluidState,
                                    int phaseIdx,
                                    int componentIdx)
    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        if (componentIdx == comp0Idx)
        {
            return MainComponent::gasEnthalpy(T, p);
        }
        else if (componentIdx == comp1Idx)
        {
            return SecondComponent::gasEnthalpy(T, p);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << componentIdx);
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the pure component at a given pressure and temperature.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar viscosity(Scalar temperature, Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "Viscosity");
    }

    using Base<Scalar, ThisType>::viscosity;
    //! \copydoc Base<Scalar,ThisType>::viscosity(const FluidState&,int)
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            const int phaseIdx)
    {
        // TODO: Add other policies especially one that allows to choose one Compoents viscosity as overall viscosity.

        // Wilke method (Reid et al.):
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        Scalar muResult = 0;
        const Scalar mu[numComponents] = {
            MainComponent::gasViscosity(T, p),
            SecondComponent::gasViscosity(T, p)
        };

        Scalar sumx = 0.0;
        using std::max;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            sumx += fluidState.moleFraction(phaseIdx, compIdx);
        sumx = max(1e-10, sumx);

        for (int i = 0; i < numComponents; ++i) {
            Scalar divisor = 0;
            using std::pow;
            using std::sqrt;
            for (int j = 0; j < numComponents; ++j) {
                Scalar phiIJ = 1 + sqrt(mu[i]/mu[j]) * pow(molarMass(j)/molarMass(i), 1/4.0);
                phiIJ *= phiIJ;
                phiIJ /= sqrt(8*(1 + molarMass(i)/molarMass(j)));
                divisor += fluidState.moleFraction(phaseIdx, j)/sumx * phiIJ;
            }
            muResult += fluidState.moleFraction(phaseIdx, i)/sumx * mu[i] / divisor;
        }
        return muResult;
    }

    using Base<Scalar, ThisType>::fugacityCoefficient;
    //! \copydoc Base<Scalar,ThisType>::fugacityCoefficient(const FluidState&,int,int)
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        if (phaseIdx == compIdx)
            // We could calculate the real fugacity coefficient of
            // the component in the fluid. Probably that's not worth
            // the effort, since the fugacity coefficient of the other
            // component is infinite anyway...
            return 1.0;
        return std::numeric_limits<Scalar>::infinity();
    }

    using Base<Scalar, ThisType>::diffusionCoefficient;
    //! \copydoc Base<Scalar,ThisType>::diffusionCoefficient(const FluidState&,int,int)
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients");
    }

    using Base<Scalar, ThisType>::binaryDiffusionCoefficient;
    //! \copydoc Base<Scalar,ThisType>::binaryDiffusionCoefficient(const FluidState&,int,int,int)
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        // TEST: set a static parameter from the input file:
        static const Scalar DiffCoeff = getParam<Scalar>("FluidSystem.DiffusionCoefficient");
        return DiffCoeff;
    }

    /*!
     * \brief Thermal conductivity of the fluid \f$\mathrm{[W/(m K)]}\f$.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar thermalConductivity(Scalar temperature, Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "Thermal Conductivity");
    }

    using Base<Scalar, ThisType>::thermalConductivity;
    //! \copydoc Base<Scalar,ThisType>::thermalConductivity(const FluidState&,int)
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "Thermal Conductivity");
    }

    /*!
     * \brief Specific isobaric heat capacity of the fluid \f$\mathrm{[J/(kg K)]}\f$.
     * \param temperature The given temperature \f$\mathrm{[K]}\f$
     * \param pressure The given pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar heatCapacity(Scalar temperature, Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "Heat Capacity");
    }

    using Base<Scalar, ThisType>::heatCapacity;
    //! \copydoc Base<Scalar,ThisType>::heatCapacity(const FluidState&,int)
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               const int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "Heat Capacity");
    }
};

} // namespace FluidSystems
} // namespace

#endif
