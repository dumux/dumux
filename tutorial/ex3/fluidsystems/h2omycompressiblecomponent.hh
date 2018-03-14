// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 *
 * \brief A fluid system with water and a fictitious component, which is to be
 *        implemented in tutorial exercise 3a, as phases and components. This
 *        fluid system is to be implemented in exercise 3b.
 */
#ifndef DUMUX_H2O_MYCOMPRESSIBLECOMPONENT_FLUID_SYSTEM_HH
#define DUMUX_H2O_MYCOMPRESSIBLECOMPONENT_FLUID_SYSTEM_HH

#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

#include <dumux/material/fluidsystems/base.hh>

// the fictitious component that was created in exercise 3a
#include <tutorial/ex3/components/mycompressiblecomponent.hh>

// the binary coefficients corresponding to this fluid system
#include <tutorial/ex3/binarycoefficients/h2omycompressiblecomponent.hh>

namespace Dumux
{
namespace FluidSystems
{

/*!
 * \brief A compositional fluid consisting of two liquid phases,
 *        which are water and a fictitious component from tutorial exercise 3a.
 */
template <class Scalar,
          class H2OType = Dumux::Components::TabulatedComponent<Dumux::Components::H2O<Scalar> > >
class H2OMyCompressibleComponent
    : public BaseFluidSystem< Scalar, H2OMyCompressibleComponent<Scalar, H2OType> >
{
    typedef H2OMyCompressibleComponent<Scalar, H2OType> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;

public:
    typedef Dumux::MyCompressibleComponent<Scalar> MyCompressibleComponent;
    typedef H2OType H2O;

    static constexpr int numPhases = 2;
    static constexpr int numComponents = 2;

    static constexpr int phase0Idx = 0; // index of the first phase
    static constexpr int phase1Idx = 1; // index of the second phase

    static constexpr int H2OIdx = 0;
    static constexpr int NAPLIdx = 1;
    // export component indices to indicate the main component
    // of the corresponding phase at atmospheric pressure 1 bar
    // and room temperature 20°C:
    static constexpr int comp0Idx = H2OIdx;
    static constexpr int comp1Idx = NAPLIdx;

    /*!
     * \brief Initialize the fluid system's static parameters generically
     *
     * If a tabulated H2O component is used, we do our best to create
     * tables that always work.
     */
    static void init()
    {
        init(/*tempMin=*/273.15,
             /*tempMax=*/623.15,
             /*numTemp=*/100,
             /*pMin=*/0.0,
             /*pMax=*/20e6,
             /*numP=*/200);
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges
     *
     * \param tempMin The minimum temperature used for tabulation of water [K]
     * \param tempMax The maximum temperature used for tabulation of water [K]
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water [Pa]
     * \param pressMax The maximum pressure used for tabulation of water [Pa]
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        if (H2O::isTabulated) {
            std::cout << "Initializing tables for the H2O fluid properties ("
                      << nTemp*nPress
                      << " entries).\n";

            H2O::init(tempMin, tempMax, nTemp,
                      pressMin, pressMax, nPress);
        }
    }


    /*!
     * \brief Return whether a phase is liquid
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isLiquid(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return true;
    }

    static constexpr bool isIdealGas(int phaseIdx)
    { return H2O::gasIsIdeal() && MyCompressibleComponent::gasIsIdeal(); }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are indepent on the fluid composition. This assumtion is true
     * if Henry's law and Raoult's law apply. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isIdealMixture(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return true;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be compressible.
     *
     * Compressible means that the partial derivative of the density
     * to the fluid pressure is always larger than zero.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isCompressible(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        if (phaseIdx == phase0Idx)
            // the water component decides for the water phase...
            return H2O::liquidIsCompressible();

        // the NAPL component decides for the napl phase...
        return MyCompressibleComponent::liquidIsCompressible();
    }

    /*!
     * \brief Return the human readable name of a phase (used in indices)
     */
    static std::string phaseName(int phaseIdx)
    {
        switch (phaseIdx) {
        case phase0Idx: return "w";
        case phase1Idx: return "n";
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return the human readable name of a component (used in indices)
     */
    static std::string componentName(int compIdx)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::name();
        case NAPLIdx: return MyCompressibleComponent::name();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the molar mass of a component in [kg/mol].
     */
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::molarMass();
        case NAPLIdx: return MyCompressibleComponent::molarMass();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Given all mole fractions in a phase, return the phase
     *        density [kg/m^3].
     */
    using Base::density;
    template <class FluidState>
    static Scalar density(const FluidState &fluidState, int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        if (phaseIdx == phase0Idx) {
            // See: doctoral thesis of Steffen Ochs 2007
            // Steam injection into saturated porous media : process analysis including experimental and numerical investigations
            // http://elib.uni-stuttgart.de/bitstream/11682/271/1/Diss_Ochs_OPUS.pdf

            // Scalar rholH2O = H2O::liquidDensity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
            // Scalar clH2O = rholH2O/H2O::molarMass();
            // Scalar x_H2O = fluidState.moleFraction (phase0Idx, H2OIdx);
            // Scalar x_myComp = fluidState.moleFraction (phase0Idx, NAPLIdx);

            /*!
             * TODO: implement the composition-dependent water density from the exercise sheet.
             */
            return /*???*/1000.0;
        }
        else {
            // assume the density of the fictious component to be independent of the composition
            Scalar pressure = MyCompressibleComponent::liquidIsCompressible()?fluidState.pressure(phaseIdx):1e100;
            return MyCompressibleComponent::liquidDensity(fluidState.temperature(phaseIdx), pressure);
        }
    }

    using Base::molarDensity;
    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density for the simple relation is defined by the
     * mass density \f$\rho_\alpha\f$ and the molar mass of the main component
     *
     * The molar density for the complrex relation is defined by the
     * mass density \f$\rho_\alpha\f$ and the mean molar mass \f$\overline M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{\overline M_\alpha} \;.\f]
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == phase0Idx)
        {
            // assume pure water or that each gas molecule displaces exactly one
            // molecule in the liquid.
            return H2O::liquidMolarDensity(T, p);
        }
        else{
            // assume the molar density of the fictious component to be independent of the composition
            Scalar pressure = MyCompressibleComponent::liquidIsCompressible()?fluidState.pressure(phaseIdx):1e100;
            return MyCompressibleComponent::liquidMolarDensity(fluidState.temperature(phaseIdx), pressure);
        }
    }

    /*!
     * \brief Return the viscosity of a phase.
     */
    using Base::viscosity;
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        if (phaseIdx == phase0Idx) {
            // assume pure water viscosity
            return H2O::liquidViscosity(fluidState.temperature(phaseIdx),
                                        fluidState.pressure(phaseIdx));
        }
        else {
            // assume pure NAPL viscosity
            return MyCompressibleComponent::liquidViscosity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        }
    }

    using Base::diffusionCoefficient;
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState, int phaseIdx, int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients");
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for components
     *        \f$\mathrm{i}\f$ and \f$\mathrm{j}\f$ in this phase.
     * \param fluidState The fluid state
     * \param paramCache mutable parameters
     * \param phaseIdx Index of the fluid phase
     * \param compIIdx Index of the component i
     * \param compJIdx Index of the component j
     */
    using Base::binaryDiffusionCoefficient;
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIIdx  && compIIdx < numComponents);
        assert(0 <= compJIdx  && compJIdx < numComponents);

        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        // we assume the diffusion coefficient to be the same in both phases
        return Dumux::BinaryCoeff::H2O_MyCompressibleComponent::liquidDiffCoeff(T, p);
    }

     /* Henry coefficients
     */
    template <class FluidState>
    static Scalar henryCoefficient(const FluidState &fluidState,
                                   int phaseIdx,
                                   int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        if (compIdx == NAPLIdx && phaseIdx == phase0Idx)
            return Dumux::BinaryCoeff::H2O_MyCompressibleComponent::henryMyCompressibleComponentInWater(T)/p;

        else if (phaseIdx == phase1Idx && compIdx == H2OIdx)
            return Dumux::BinaryCoeff::H2O_MyCompressibleComponent::henryWaterInMyCompressibleComponent(T)/p;

        else
            DUNE_THROW(Dune::InvalidStateException, "non-existent henry coefficient for phase index " << phaseIdx
                                                     << " and component index " << compIdx);
    }

    using Base::fugacityCoefficient;
    /*!
     * \brief Returns the fugacity coefficient \f$\mathrm{[-]}\f$ of a component in a
     *        phase.
     *
     * In this case, things are actually pretty simple. We have an ideal
     * solution. Thus, the fugacity coefficient is 1 in the gas phase
     * (fugacity equals the partial pressure of the component in the gas phase
     * respectively in the liquid phases it is the inverse of the
     * Henry coefficients scaled by pressure
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase
     * \param compIdx The index of the component
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == phase0Idx) {
            if (compIdx == H2OIdx)
                return H2O::vaporPressure(T)/p;
            else if (compIdx == NAPLIdx)
                return Dumux::BinaryCoeff::H2O_MyCompressibleComponent::henryMyCompressibleComponentInWater(T)/p;
        }

        // for the NAPL phase, we assume currently that nothing is
        // dissolved. this means that the affinity of the NAPL
        // component to the NAPL phase is much higher than for the
        // other components, i.e. the fugacity coefficient is much
        // smaller.
        Scalar phiNapl = MyCompressibleComponent::vaporPressure(T)/p;
        if (compIdx == NAPLIdx)
            return phiNapl;
        else
            return 1e6*phiNapl;
    }

    template <class FluidState>
    static Scalar kelvinVaporPressure(const FluidState &fluidState,
                                      const int phaseIdx,
                                      const int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::H2OMyCompressibleComponent::kelvinVaporPressure()");
    }

     /*  partial pressures in the gas phase, taken from saturation vapor pressures
     */
    template <class FluidState>
    static Scalar partialPressureGas(const FluidState &fluidState, int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= compIdx  && compIdx < numComponents);

        const Scalar T = fluidState.temperature(phaseIdx);
        if (compIdx == NAPLIdx)
            return MyCompressibleComponent::vaporPressure(T);
        else if (compIdx == H2OIdx)
            return H2O::vaporPressure(T);
        else
            DUNE_THROW(Dune::InvalidStateException, "non-existent component index " << compIdx);
    }

     /*  inverse vapor pressures, taken from inverse saturation vapor pressures
     */
    template <class FluidState>
    static Scalar inverseVaporPressureCurve(const FluidState &fluidState,
                                            int phaseIdx,
                                            int compIdx)
    {
        assert(0 <= compIdx  && compIdx < numComponents);

        const Scalar pressure = fluidState.pressure(phaseIdx);
        if (compIdx == NAPLIdx)
            return MyCompressibleComponent::vaporTemperature(pressure);
        else if (compIdx == H2OIdx)
            return H2O::vaporTemperature(pressure);
        else
            DUNE_THROW(Dune::InvalidStateException, "non-existent component index " << compIdx);
    }



    /*!
     * \brief Given all mole fractions in a phase, return the specific
     *        phase enthalpy [J/kg].
     */
    using Base::enthalpy;
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::H2OMyCompressibleComponent::enthalpy()");
    }

    using Base::heatCapacity;
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::H2OMyCompressibleComponent::heatCapacity()");
    }

    using Base::thermalConductivity;
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const Scalar temperature  = fluidState.temperature(phaseIdx) ;
        const Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == phase0Idx)
        {
            return H2O::liquidThermalConductivity(temperature, pressure);
        }
        else
        {
            return MyCompressibleComponent::liquidThermalConductivity(temperature, pressure);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

private:

};
} // end namespace FluidSystems
} // end namespace Dumux

#endif
