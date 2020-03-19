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
 * \brief @copybrief Dumux::FluidSystems::H2OHeavyOil
 */
#ifndef DUMUX_H2O_HEAVYOIL_FLUID_SYSTEM_HH
#define DUMUX_H2O_HEAVYOIL_FLUID_SYSTEM_HH

#include <dumux/material/idealgas.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/heavyoil.hh>

#include <dumux/material/binarycoefficients/h2o_heavyoil.hh>

#include <dumux/material/fluidsystems/base.hh>

#include <dumux/io/name.hh>

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 * \brief A compositional fluid system with water and heavy oil
 *        components in both the liquid and the gas phase.
 */
template <class Scalar,
          class H2OType = Dumux::Components::TabulatedComponent<Dumux::Components::H2O<Scalar> > >
class H2OHeavyOil
    : public Base<Scalar, H2OHeavyOil<Scalar, H2OType> >
{
    using ThisType = H2OHeavyOil<Scalar, H2OType>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;

public:
    using HeavyOil = Dumux::Components::HeavyOil<Scalar>;
    using H2O = H2OType;


    static const int numPhases = 3;
    static const int numComponents = 2;

    static const int wPhaseIdx = 0; // index of the water phase
    static const int nPhaseIdx = 1; // index of the NAPL phase
    static const int gPhaseIdx = 2; // index of the gas phase

    static const int H2OIdx = 0;
    static const int NAPLIdx = 1;

    // export component indices to indicate the main component
    // of the corresponding phase at atmospheric pressure 1 bar
    // and room temperature 20°C:
    static const int wCompIdx = H2OIdx;
    static const int nCompIdx = NAPLIdx;

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
        if (H2O::isTabulated)
        {
            H2O::init(tempMin, tempMax, nTemp,
                      pressMin, pressMax, nPress);
        }
    }

    /*!
     * \brief Get the main component of a given phase
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr int getMainComponent(int phaseIdx)
    {
        // For the gas phase, choosing a main component appears to be
        // rather arbitrary. Motivated by the fact that the thermal conductivity
        // of the gas phase is set to the thermal conductivity of pure water,
        // water is chosen for now.
        if (phaseIdx == nPhaseIdx)
            return nCompIdx;
        else
            return wCompIdx;
    }

    /*!
     * \brief Returns whether the fluids are miscible
     */
    static constexpr bool isMiscible()
    { return true; }

    /*!
     * \brief Return whether a phase is gaseous
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx == gPhaseIdx;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal gas.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isIdealGas(int phaseIdx)
    { return phaseIdx == gPhaseIdx && H2O::gasIsIdeal() && HeavyOil::gasIsIdeal(); }

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
    static bool isIdealMixture(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        // we assume Henry's and Raoult's laws for the water phase and
        // and no interaction between gas molecules of different
        // components, so all phases are ideal mixtures!
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
        // gases are always compressible
        if (phaseIdx == gPhaseIdx)
            return true;
        else if (phaseIdx == wPhaseIdx)
            // the water component decides for the water phase...
            return H2O::liquidIsCompressible();

        // the NAPL component decides for the napl phase...
        return HeavyOil::liquidIsCompressible();
    }

    /*!
     * \brief Return the human readable name of a phase (used in indices)
     */
    static std::string phaseName(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        switch (phaseIdx)
        {
            case wPhaseIdx: return IOName::aqueousPhase();
            case nPhaseIdx: return IOName::naplPhase();
            case gPhaseIdx: return IOName::gaseousPhase();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return the human readable name of a component (used in indices)
     */
    static std::string componentName(int compIdx)
    {
        switch (compIdx) {
            case H2OIdx: return H2O::name();
            case NAPLIdx: return HeavyOil::name();
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
            case NAPLIdx: return HeavyOil::molarMass();
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
        if (phaseIdx == wPhaseIdx) {
            // See: doctoral thesis of Steffen Ochs 2007
            // Steam injection into saturated porous media : process analysis including experimental and numerical investigations
            // http://elib.uni-stuttgart.de/bitstream/11682/271/1/Diss_Ochs_OPUS.pdf

            // This assumes each gas molecule displaces exactly one
            // molecule in the liquid.

            return H2O::liquidMolarDensity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx))
                   * (H2O::molarMass()*fluidState.moleFraction(wPhaseIdx, H2OIdx)
                      + HeavyOil::molarMass()*fluidState.moleFraction(wPhaseIdx, NAPLIdx));
        }
        else if (phaseIdx == nPhaseIdx) {
            // assume pure NAPL for the NAPL phase
            Scalar pressure = HeavyOil::liquidIsCompressible()?fluidState.pressure(phaseIdx):1e100;
            return HeavyOil::liquidDensity(fluidState.temperature(phaseIdx), pressure);
        }

        assert (phaseIdx == gPhaseIdx);
        Scalar pH2O =
            fluidState.moleFraction(gPhaseIdx, H2OIdx)  *
            fluidState.pressure(gPhaseIdx);
        Scalar pNAPL =
            fluidState.moleFraction(gPhaseIdx, NAPLIdx)  *
            fluidState.pressure(gPhaseIdx);
        return H2O::gasDensity(fluidState.temperature(phaseIdx), pH2O)
               + HeavyOil::gasDensity(fluidState.temperature(phaseIdx), pNAPL);
    }

    using Base::molarDensity;
    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density is defined by the
     * mass density \f$\rho_\alpha\f$ and the mean molar mass \f$\overline M_\alpha\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{\overline M_\alpha} \;.\f]
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, int phaseIdx)
    {
        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == nPhaseIdx)
        {
            return HeavyOil::liquidMolarDensity(temperature, pressure);
        }
        else if (phaseIdx == wPhaseIdx)
        {   // This assumes each gas molecule displaces exactly one
            // molecule in the liquid.
            return H2O::liquidMolarDensity(temperature, pressure);
        }
        else
        {
            return H2O::gasMolarDensity(temperature, fluidState.partialPressure(gPhaseIdx, H2OIdx))
                   + HeavyOil::gasMolarDensity(temperature, fluidState.partialPressure(gPhaseIdx, NAPLIdx));
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
        if (phaseIdx == wPhaseIdx) {
            // assume pure water viscosity
            return H2O::liquidViscosity(fluidState.temperature(phaseIdx),
                                        fluidState.pressure(phaseIdx));
        }
        else if (phaseIdx == nPhaseIdx) {
            // assume pure NAPL viscosity
            return HeavyOil::liquidViscosity(fluidState.temperature(phaseIdx),
                                             fluidState.pressure(phaseIdx));
        }

        assert (phaseIdx == gPhaseIdx);

        /* Wilke method. See:
         *
         * See: R. Reid, et al.: The Properties of Gases and Liquids,
         * 4th edition, McGraw-Hill, 1987, 407-410
         * 5th edition, McGraw-Hill, 20001, p. 9.21/22
         *
         * in this case, we use a simplified version in order to avoid
         * computationally costly evaluation of sqrt and pow functions and
         * divisions
         * -- compare e.g. with Promo Class p. 32/33
         */
        const Scalar mu[numComponents] = {
            H2O::gasViscosity(fluidState.temperature(phaseIdx), H2O::vaporPressure(fluidState.temperature(phaseIdx))),
            HeavyOil::gasViscosity(fluidState.temperature(phaseIdx), HeavyOil::vaporPressure(fluidState.temperature(phaseIdx)))
        };

        return mu[H2OIdx]*fluidState.moleFraction(gPhaseIdx, H2OIdx)
               + mu[NAPLIdx]*fluidState.moleFraction(gPhaseIdx, NAPLIdx);
    }

    using Base::binaryDiffusionCoefficient;
    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for components
     *        \f$\mathrm{i}\f$ and \f$\mathrm{j}\f$ in this phase.
     *
     * Be aware that in this case there are only two components in three phases.
     * Therefore, we assume the diffusion to simply be the binary diffusion coefficients.
     * This was previously implemented in diffusionCoefficient(), but is now moved.
     *
     * \param fluidState The fluid state
     * \param phaseIdx Index of the fluid phase
     * \param compIIdx Index of the component i
     * \param compJIdx Index of the component j
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)

    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == wPhaseIdx)
            return BinaryCoeff::H2O_HeavyOil::liquidDiffCoeff(T, p);

        // gas phase
        else if (phaseIdx == gPhaseIdx)
            return BinaryCoeff::H2O_HeavyOil::gasDiffCoeff(T, p);
        else
            DUNE_THROW(Dune::InvalidStateException,
                       "Non-existent binary diffusion coefficient for phase index "
                       << phaseIdx);

    }

    using Base::diffusionCoefficient;
    /*!
     * \brief Calculate the binary molecular diffusion coefficient for
     *        a component in a fluid phase \f$\mathrm{[mol^2 * s / (kg*m^3)]}\f$
     * \param fluidState The fluid state
     * \param phaseIdx Index of the fluid phase
     * Molecular diffusion of a component \f$\mathrm{\kappa}\f$ is caused by a
     * gradient of the chemical potential and follows the law
     *
     * \f[ J = - D \mathbf{grad} \mu_\kappa \f]
     *
     * where \f$\mathrm{\mu_\kappa}\f$ is the component's chemical potential,
     * \f$\mathrm{D}\f$ is the diffusion coefficient and \f$\mathrm{J}\f$ is the
     * diffusive flux. \f$\mathrm{\mu_\kappa}\f$ is connected to the component's
     * fugacity \f$\mathrm{f_\kappa}\f$ by the relation
     *
     * \f[ \mu_\kappa = R T_\alpha \mathrm{ln} \frac{f_\kappa}{p_\alpha} \f]
     *
     * where \f$\mathrm{p_\alpha}\f$ and \f$\mathrm{T_\alpha}\f$ are the fluid phase'
     * pressure and temperature.
     */
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState, int phaseIdx)
    {
        // liquid phase
        if (phaseIdx == wPhaseIdx)
            return binaryDiffusionCoefficient(fluidState, phaseIdx, H2OIdx, NAPLIdx);
        // gas phase
        else if (phaseIdx == gPhaseIdx)
            return binaryDiffusionCoefficient(fluidState, phaseIdx, NAPLIdx, H2OIdx);
        else
            DUNE_THROW(Dune::InvalidStateException,
                       "Non-existent diffusion coefficient for phase index "<< phaseIdx);
    }

    /*!
     * \brief Henry coefficients \f$[N/m^2]\f$ of a component in a phase.
     */
    template <class FluidState>
    static Scalar henryCoefficient(const FluidState &fluidState,
                                   int phaseIdx,
                                   int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        const Scalar T = fluidState.temperature(phaseIdx);

        if (compIdx == NAPLIdx && phaseIdx == wPhaseIdx)
            return Dumux::BinaryCoeff::H2O_HeavyOil::henryOilInWater(T);

        else if (phaseIdx == nPhaseIdx && compIdx == H2OIdx)
            return Dumux::BinaryCoeff::H2O_HeavyOil::henryWaterInOil(T);

        else
            DUNE_THROW(Dune::InvalidStateException, "non-existent henry coefficient for phase index " << phaseIdx
                                                     << " and component index " << compIdx);
    }

    /*!
     * \brief Partial pressures in the gas phase, taken from saturation vapor pressures.
     */
    template <class FluidState>
    static Scalar partialPressureGas(const FluidState &fluidState, int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= compIdx  && compIdx < numComponents);

        const Scalar T = fluidState.temperature(phaseIdx);
        if (compIdx == NAPLIdx)
            return HeavyOil::vaporPressure(T);
        else if (compIdx == H2OIdx)
            return H2O::vaporPressure(T);
        else
            DUNE_THROW(Dune::InvalidStateException, "non-existent component index " << compIdx);
    }

    /*!
     * \brief Inverse vapor pressures, taken from inverse saturation vapor pressures
     */
    template <class FluidState>
    static Scalar inverseVaporPressureCurve(const FluidState &fluidState,
                                            int phaseIdx,
                                            int compIdx)
    {
        assert(0 <= compIdx  && compIdx < numComponents);

        const Scalar pressure = fluidState.pressure(phaseIdx);
        if (compIdx == NAPLIdx)
            return HeavyOil::vaporTemperature(pressure);
        else if (compIdx == H2OIdx)
            return H2O::vaporTemperature(pressure);
        else
            DUNE_THROW(Dune::InvalidStateException, "non-existent component index " << compIdx);
    }



    /*!
     * \brief Given all mole fractions in a phase, return the specific
     *        phase enthalpy\f$\mathrm{[J/kg]}\f$.
     *  \todo This system neglects the contribution of gas-molecules in the liquid phase.
     *        This contribution is probably not big.
     *        Somebody would have to find out the enthalpy of solution for this system. ...
     */
    using Base::enthalpy;
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           int phaseIdx)
    {
        if (phaseIdx == wPhaseIdx) {
            return H2O::liquidEnthalpy(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        }
        else if (phaseIdx == nPhaseIdx) {
            return HeavyOil::liquidEnthalpy(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        }
        else if (phaseIdx == gPhaseIdx) {  // gas phase enthalpy depends strongly on composition
            Scalar hgc = HeavyOil::gasEnthalpy(fluidState.temperature(phaseIdx),
                                           fluidState.pressure(phaseIdx));
            Scalar hgw = H2O::gasEnthalpy(fluidState.temperature(phaseIdx),
                                          fluidState.pressure(phaseIdx));

            Scalar result = 0;
            result += hgw * fluidState.massFraction(gPhaseIdx, H2OIdx);
            result += hgc * fluidState.massFraction(gPhaseIdx, NAPLIdx);

            return result;
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Returns the specific enthalpy \f$\mathrm{[J/kg]}\f$ of a component in a specific phase
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase
     * \param componentIdx The index of the component
     */
    template <class FluidState>
    static Scalar componentEnthalpy(const FluidState& fluidState, int phaseIdx, int componentIdx)
    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == wPhaseIdx)
        {
            if (componentIdx == H2OIdx)
                return H2O::liquidEnthalpy(T, p);
            else if (componentIdx == NAPLIdx)
                DUNE_THROW(Dune::NotImplemented, "The component enthalpy for NAPL in water is not implemented.");
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << componentIdx);
        }
        else if (phaseIdx == nPhaseIdx)
        {
            if (componentIdx == H2OIdx)
                DUNE_THROW(Dune::NotImplemented, "The component enthalpy for water in NAPL is not implemented.");
            else if (componentIdx == NAPLIdx)
                return HeavyOil::liquidEnthalpy(T, p);
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << componentIdx);
        }
        else if (phaseIdx == gPhaseIdx)
        {
            if (componentIdx == H2OIdx)
                return H2O::gasEnthalpy(T, p);
            else if (componentIdx == NAPLIdx)
                return HeavyOil::gasEnthalpy(T, p);
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << componentIdx);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    using Base::heatCapacity;
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::H2ONAPL::heatCapacity()");
    }

    /*!
     * \brief Thermal conductivity of a fluid phase \f$\mathrm{[W/(m K)]}\f$.
     *
     * Use the conductivity of water (wPhase and gPhase) and oil (nPhase) as a first approximation.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    using Base::thermalConductivity;
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      int phaseIdx)
    {
        const Scalar temperature  = fluidState.temperature(phaseIdx) ;
        const Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == wPhaseIdx)
        {
            return H2O::liquidThermalConductivity(temperature, pressure);
        }
        else if (phaseIdx == nPhaseIdx)
        {
            return HeavyOil::liquidThermalConductivity(temperature, pressure);
        }
        else if (phaseIdx == gPhaseIdx)
        {
            return H2O::gasThermalConductivity(temperature, pressure);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

};
} // end namespace FluidSystems
} // end namespace Dumux

#endif
