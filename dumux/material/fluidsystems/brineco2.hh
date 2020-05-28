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
 * \brief @copybrief Dumux::FluidSystems::BrineCO2
 */
#ifndef DUMUX_BRINE_CO2_FLUID_SYSTEM_HH
#define DUMUX_BRINE_CO2_FLUID_SYSTEM_HH

#include <type_traits>

#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>
#include <dumux/material/idealgas.hh>
#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/fluidsystems/brine.hh>
#include <dumux/material/fluidstates/adapter.hh>

#include <dumux/material/components/brine.hh>
#include <dumux/material/components/co2.hh>
#include <dumux/material/components/co2tablereader.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

#include <dumux/material/binarycoefficients/brine_co2.hh>

#include <dumux/io/name.hh>

namespace Dumux {

// include the default tables for CO2
#ifndef DOXYGEN // hide tables from doxygen
#include <dumux/material/components/co2tables.inc>
#endif

namespace FluidSystems {
namespace Detail {

    /*!
     * \brief Class that exports some indices that should
     *        be provided by the BrineAir fluid system.
     *        The indices are chosen dependent on the policy,
     *        i.e. if a simplified pseudo component Brine is
     *        used or salt is considered an individual component.
     */
    template<bool useConstantSalinity>
    struct BrineCO2Indices;

    /*!
     * \brief Specialization for the case of brine being
     *        a pseudo component with a constant salinity.
     * \note This specialization exports brine as component
     */
    template<>
    struct BrineCO2Indices<true>
    {
        static constexpr int BrineIdx = 0;
    };

    /*!
     * \brief Specialization for the case of brine being
     *        a fluid system with NaCl as individual component.
     * \note This specialization exports both water and NaCl as components
     */
    template<>
    struct BrineCO2Indices<false>
    {
        static constexpr int H2OIdx = 0;
        static constexpr int NaClIdx = 2;
        static constexpr int comp2Idx = 2;
    };
} // end namespace Detail

/*!
 * \ingroup Fluidsystems
 * \brief Default policy for the Brine-CO2 fluid system
 */
template<bool salinityIsConstant, bool fastButSimplifiedRelations = false>
struct BrineCO2DefaultPolicy
{
    static constexpr bool useConstantSalinity() { return salinityIsConstant; }
    static constexpr bool useCO2GasDensityAsGasMixtureDensity() { return fastButSimplifiedRelations; }
};

/*!
 * \ingroup Fluidsystems
 * \brief A compositional fluid with brine (H2O & NaCl) and carbon dioxide as
 *        components in both the liquid and the gas (supercritical) phase.
 *
 * \note Depending on the chosen policy, the salinity is assumed to be constant
 *       (in which case Brine is used as a pseudo component) or salt (here NaCl)
 *       is considered as an individual component.
 * \note This implemetation always assumes NaCl stays in the liquid phase.
 */
template< class Scalar,
          class CO2Table,
          class H2OType = Components::TabulatedComponent<Components::H2O<Scalar>>,
          class Policy = BrineCO2DefaultPolicy</*constantSalinity?*/true> >
class BrineCO2
: public Base<Scalar, BrineCO2<Scalar, CO2Table, H2OType, Policy>>
, public Detail::BrineCO2Indices<Policy::useConstantSalinity()>
{
    using ThisType = BrineCO2<Scalar, CO2Table, H2OType, Policy>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;

    // binary coefficients
    using Brine_CO2 = BinaryCoeff::Brine_CO2<Scalar, CO2Table>;

    // use constant salinity brine?
    static constexpr bool useConstantSalinity = Policy::useConstantSalinity();

    // The possible brine types
    using VariableSalinityBrine = Dumux::FluidSystems::Brine<Scalar, H2OType>;
    using ConstantSalinityBrine = Dumux::Components::Brine<Scalar, H2OType>;
    using BrineType = typename std::conditional_t< useConstantSalinity,
                                                   ConstantSalinityBrine,
                                                   VariableSalinityBrine >;

    /////////////////////////////////////////////////////////////////////////////////
    //! The following two indices are only used internally and are not part of the
    //! public interface. Depending on the chosen policy, i.e. if brine is used as
    //! a pseudo component or a fluid system with NaCl as a separate component, the
    //! indices that are part of the public interface are chosen by inheritance from
    //! Detail::BrineCO2Indices (see documentation).
    //!
    //! depending on the implementation this is either brine (pseudo-component) or H2O
    static constexpr int BrineOrH2OIdx = 0;
    //! if the implementation considers NaCl as a real compoent, it gets the index 2
    static constexpr int NaClIdx = 2;

public:
    using ParameterCache = NullParameterCache;

    using H2O = H2OType;
    using Brine = BrineType;
    using CO2 = Dumux::Components::CO2<Scalar, CO2Table>;

    static constexpr int numComponents = useConstantSalinity ? 2 : 3;
    static constexpr int numPhases = 2;

    static constexpr int liquidPhaseIdx = 0; //!< index of the liquid phase
    static constexpr int gasPhaseIdx = 1;    //!< index of the gas phase
    static constexpr int phase0Idx = liquidPhaseIdx; //!< index of the first phase
    static constexpr int phase1Idx = gasPhaseIdx;    //!< index of the second phase

    static constexpr int comp0Idx = 0;
    static constexpr int comp1Idx = 1;

    // CO2 is always the second component
    static constexpr int CO2Idx = comp1Idx;

private:

    // Adapter policy for the fluid state corresponding to the brine fluid system
    struct BrineAdapterPolicy
    {
        using FluidSystem = VariableSalinityBrine;

        static constexpr int phaseIdx(int brinePhaseIdx) { return liquidPhaseIdx; }
        static constexpr int compIdx(int brineCompIdx)
        {
            switch (brineCompIdx)
            {
                assert(brineCompIdx == VariableSalinityBrine::H2OIdx || brineCompIdx == VariableSalinityBrine::NaClIdx);
                case VariableSalinityBrine::H2OIdx: return BrineOrH2OIdx;
                case VariableSalinityBrine::NaClIdx: return NaClIdx;
                default: return 0; // this will never be reached, only needed to suppress compiler warning
            }
        }
    };

    template<class FluidState>
    using BrineAdapter = FluidStateAdapter<FluidState, BrineAdapterPolicy>;

public:

    /*!
     * \brief Return the human readable name of a fluid phase
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx)
    {
        switch (phaseIdx)
        {
            case liquidPhaseIdx: return IOName::liquidPhase();
            case gasPhaseIdx: return IOName::gaseousPhase();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
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

        return phaseIdx == gasPhaseIdx;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal gas.
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isIdealGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        // let the fluids decide
        if (phaseIdx == gasPhaseIdx)
            return useConstantSalinity ? (ConstantSalinityBrine::gasIsIdeal() && CO2::gasIsIdeal())
                                       : (H2O::gasIsIdeal() && CO2::gasIsIdeal());
        return false; // not a gas
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are independent on the fluid composition. This assumption is true
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
        if (phaseIdx == liquidPhaseIdx)
            return false;
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
        if (phaseIdx == liquidPhaseIdx)
            return useConstantSalinity ? ConstantSalinityBrine::liquidIsCompressible()
                                       : VariableSalinityBrine::isCompressible(VariableSalinityBrine::liquidPhaseIdx);
        return true;
    }

    /*!
     * \brief Return the human readable name of a component
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);
        if (useConstantSalinity)
        {
            static std::string name[] = { ConstantSalinityBrine::name(), CO2::name() };
            return name[compIdx];
        }
        else
        {
            static std::string name[] = { VariableSalinityBrine::componentName(VariableSalinityBrine::H2OIdx),
                                          CO2::name(),
                                          VariableSalinityBrine::NaCl::name() };
            return name[compIdx];
        }
    }

    /*!
     * \brief Return the molar mass of a component in \f$\mathrm{[kg/mol]}\f$.
     * \param compIdx The index of the component to consider
     */
    static Scalar molarMass(int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);
        if (useConstantSalinity)
        {
            static const Scalar M[] = { ConstantSalinityBrine::molarMass(), CO2::molarMass() };
            return M[compIdx];
        }
        else
        {
            static const Scalar M[] = { VariableSalinityBrine::molarMass(VariableSalinityBrine::H2OIdx),
                                        CO2::molarMass(),
                                        VariableSalinityBrine::molarMass(VariableSalinityBrine::NaClIdx) };
            return M[compIdx];
        }
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    // Initializing with salinity and default tables
    static void init()
    {
        init(/*startTemp=*/273.15, /*endTemp=*/623.15, /*tempSteps=*/100,
             /*startPressure=*/1e4, /*endPressure=*/40e6, /*pressureSteps=*/200);
    }

    // Initializing with custom tables
    static void init(Scalar startTemp, Scalar endTemp, int tempSteps,
                     Scalar startPressure, Scalar endPressure, int pressureSteps)
    {
        std::cout << "The Brine-CO2 fluid system was configured with the following policy:\n";
        std::cout << " - use constant salinity: " << std::boolalpha << Policy::useConstantSalinity() << "\n";
        std::cout << " - use CO2 gas density as gas mixture density: " << std::boolalpha << Policy::useCO2GasDensityAsGasMixtureDensity() << std::endl;

        if (H2O::isTabulated)
            H2O::init(startTemp, endTemp, tempSteps, startPressure, endPressure, pressureSteps);
    }

    using Base::density;
    /*!
     * \brief Given a phase's composition, temperature, pressure, and
     *        the partial pressures of all components, return its
     *        density \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase
     */
    template <class FluidState>
    static Scalar density(const FluidState& fluidState, int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        if (phaseIdx == liquidPhaseIdx)
            return liquidDensityMixture_(fluidState);

        else if (phaseIdx == gasPhaseIdx)
        {
            if (Policy::useCO2GasDensityAsGasMixtureDensity())
                // use the CO2 gas density only and neglect compositional effects
                return CO2::gasDensity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
            else
            {
                // assume ideal mixture: steam and CO2 don't "see" each other
                Scalar rho_gH2O = H2O::gasDensity(T, fluidState.partialPressure(gasPhaseIdx, BrineOrH2OIdx));
                Scalar rho_gCO2 = CO2::gasDensity(T, fluidState.partialPressure(gasPhaseIdx, CO2Idx));
                return (rho_gH2O + rho_gCO2);
            }
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index.");
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
    static Scalar molarDensity(const FluidState& fluidState, int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        if (phaseIdx == liquidPhaseIdx)
            return density(fluidState, phaseIdx)/fluidState.averageMolarMass(phaseIdx);
        else if (phaseIdx == gasPhaseIdx)
        {
            if (Policy::useCO2GasDensityAsGasMixtureDensity())
                return CO2::gasMolarDensity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
            else
            {
                // assume ideal mixture: steam and CO2 don't "see" each other
                Scalar rhoMolar_gH2O = H2O::gasMolarDensity(T, fluidState.partialPressure(gasPhaseIdx, BrineOrH2OIdx));
                Scalar rhoMolar_gCO2 = CO2::gasMolarDensity(T, fluidState.partialPressure(gasPhaseIdx, CO2Idx));
                return rhoMolar_gH2O + rhoMolar_gCO2;
            }
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index.");
    }

    using Base::viscosity;
    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase \f$\mathrm{[Pa*s]}\f$
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     *
     * \note For the viscosity of the phases the contribution of the minor
     *       component is neglected. This contribution is probably not big, but somebody
     *       would have to find out its influence.
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState& fluidState, int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == liquidPhaseIdx)
            return useConstantSalinity ? ConstantSalinityBrine::liquidViscosity(T, p)
                                       : VariableSalinityBrine::viscosity( BrineAdapter<FluidState>(fluidState),
                                                                           VariableSalinityBrine::liquidPhaseIdx );
        else if (phaseIdx == gasPhaseIdx)
            return CO2::gasViscosity(T, p);

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index.");
    }

    using Base::fugacityCoefficient;
    /*!
     * \brief Returns the fugacity coefficient \f$\mathrm{[-]}\f$ of a component in a
     *        phase.
     *
     * The fugacity coefficient \f$\mathrm{\phi^\kappa_\alpha}\f$ of
     * component \f$\mathrm{\kappa}\f$ in phase \f$\mathrm{\alpha}\f$ is connected to
     * the fugacity \f$\mathrm{f^\kappa_\alpha}\f$ and the component's mole
     * fraction \f$\mathrm{x^\kappa_\alpha}\f$ by means of the relation
     *
     * \f[
     f^\kappa_\alpha = \phi^\kappa_\alpha\;x^\kappa_\alpha\;p_\alpha
     \f]
     * where \f$\mathrm{p_\alpha}\f$ is the pressure of the fluid phase.
     *
     * The fugacity itself is just an other way to express the
     * chemical potential \f$\mathrm{\zeta^\kappa_\alpha}\f$ of the component:
     *
     * \f[
     f^\kappa_\alpha := \exp\left\{\frac{\zeta^\kappa_\alpha}{k_B T_\alpha} \right\}
     \f]
     * where \f$\mathrm{k_B}\f$ is Boltzmann's constant.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState& fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= compIdx && compIdx < numComponents);

        if (phaseIdx == gasPhaseIdx)
            // use the fugacity coefficients of an ideal gas. the
            // actual value of the fugacity is not relevant, as long
            // as the relative fluid compositions are observed,
            return 1.0;

        else if (phaseIdx == liquidPhaseIdx)
        {
            Scalar T = fluidState.temperature(phaseIdx);
            Scalar pl = fluidState.pressure(liquidPhaseIdx);
            Scalar pg = fluidState.pressure(gasPhaseIdx);

            assert(T > 0);
            assert(pl > 0 && pg > 0);

            // calulate the equilibrium composition for given T & p
            Scalar xlH2O, xgH2O;
            Scalar xlCO2, xgCO2;
            const Scalar salinity = useConstantSalinity ? ConstantSalinityBrine::salinity()
                                                        : fluidState.massFraction(liquidPhaseIdx, NaClIdx);
            Brine_CO2::calculateMoleFractions(T, pl, salinity, /*knownGasPhaseIdx=*/-1, xlCO2, xgH2O);

            // normalize the phase compositions
            using std::min;
            using std::max;
            xlCO2 = max(0.0, min(1.0, xlCO2));
            xgH2O = max(0.0, min(1.0, xgH2O));
            xlH2O = 1.0 - xlCO2;
            xgCO2 = 1.0 - xgH2O;

            if (compIdx == BrineOrH2OIdx)
                return (xgH2O/xlH2O)*(pg/pl);

            else if (compIdx == CO2Idx)
                return (xgCO2/xlCO2)*(pg/pl);

            // NaCl is assumed to stay in the liquid!
            else if (!useConstantSalinity && compIdx == NaClIdx)
                return 0.0;

            DUNE_THROW(Dune::InvalidStateException, "Invalid component index.");
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index.");
    }

    /*!
     * \brief Returns the equilibrium concentration of the dissolved component
     *        in a phase.
     * \param fluidState An arbitrary fluid state
     * \param paramCache Parameter cache
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar equilibriumMoleFraction(const FluidState& fluidState,
                                          const ParameterCache& paramCache,
                                          int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        assert(T > 0);
        assert(p > 0);

        Scalar xgH2O;
        Scalar xlCO2;

        // calulate the equilibrium composition for given T & p
        const Scalar salinity = useConstantSalinity ? ConstantSalinityBrine::salinity()
                                                    : fluidState.massFraction(liquidPhaseIdx, NaClIdx);
        Brine_CO2::calculateMoleFractions(T, p, salinity, /*knowgasPhaseIdx=*/-1, xlCO2, xgH2O);

        if (phaseIdx == gasPhaseIdx)
            return xgH2O;
        else if (phaseIdx == liquidPhaseIdx)
            return xlCO2;

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index.");
    }


    using Base::diffusionCoefficient;
    /*!
     * \brief Calculate the molecular diffusion coefficient for a
     *        component in a fluid phase \f$\mathrm{[mol^2 * s / (kg*m^3)]}\f$
     *
     * Molecular diffusion of a compoent \f$\mathrm{\kappa}\f$ is caused by a
     * gradient of the chemical potential and follows the law
     *
     * \f[ J = - D \textbf{grad} mu_\kappa \f]
     *
     * where \f$\mathrm{\mu_\kappa}\f$ is the component's chemical potential,
     * \f$D\f$ is the diffusion coefficient and \f$\mathrm{J}\f$ is the
     * diffusive flux. \f$\mathrm{mu_\kappa}\f$ is connected to the component's
     * fugacity \f$\mathrm{f_\kappa}\f$ by the relation
     *
     * \f[ \mu_\kappa = R T_\alpha \mathrm{ln} \frac{f_\kappa}{p_\alpha} \f]
     *
     * where \f$\mathrm{p_\alpha}\f$ and \f$\mathrm{T_\alpha}\f$ are the fluid phase'
     * pressure and temperature.
     *
     * Maybe see http://www.ddbst.de/en/EED/PCP/DIF_C1050.php
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState& fluidState, int phaseIdx, int compIdx)
    { DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients"); }

    using Base::binaryDiffusionCoefficient;
    /*!
     * \brief Given the phase compositions, return the binary
     *        diffusion coefficient \f$\mathrm{[m^2/s]}\f$ of two components in a phase.
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIIdx Index of the component i
     * \param compJIdx Index of the component j
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState& fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        assert(0 <= compIIdx && compIIdx < numComponents);
        assert(0 <= compJIdx && compJIdx < numComponents);

        if (compIIdx > compJIdx)
        {
            using std::swap;
            swap(compIIdx, compJIdx);
        }

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);
        if (phaseIdx == liquidPhaseIdx)
        {
            if (compIIdx == BrineOrH2OIdx && compJIdx == CO2Idx)
                return Brine_CO2::liquidDiffCoeff(T, p);
            if (!useConstantSalinity && compIIdx == BrineOrH2OIdx && compJIdx == NaClIdx)
                return VariableSalinityBrine::binaryDiffusionCoefficient( BrineAdapter<FluidState>(fluidState),
                                                                          VariableSalinityBrine::liquidPhaseIdx,
                                                                          VariableSalinityBrine::H2OIdx,
                                                                          VariableSalinityBrine::NaClIdx );

            DUNE_THROW(Dune::NotImplemented, "Binary diffusion coefficient of components " <<
                                             compIIdx << " and " << compJIdx  << " in phase " << phaseIdx);
        }
        else if (phaseIdx == gasPhaseIdx)
        {
            if (compIIdx == BrineOrH2OIdx && compJIdx == CO2Idx)
                return Brine_CO2::gasDiffCoeff(T, p);

            // NaCl is expected to never be present in the gas phase. we need to
            // return a diffusion coefficient that does not case numerical problems.
            // We choose a very small value here.
            else if (!useConstantSalinity && compIIdx == CO2Idx && compJIdx == NaClIdx)
                return 1e-12;

            DUNE_THROW(Dune::NotImplemented, "Binary diffusion coefficient of components " <<
                                             compIIdx << " and " << compJIdx  << " in phase " << phaseIdx);
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index.");
    }

    using Base::enthalpy;
    /*!
     * \brief Given the phase composition, return the specific
     *        phase enthalpy \f$\mathrm{[J/kg]}\f$.
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState& fluidState, int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == liquidPhaseIdx)
        {
            // Convert J/kg to kJ/kg
            const Scalar h_ls1 = useConstantSalinity ? ConstantSalinityBrine::liquidEnthalpy(T, p)/1e3
                                                     : VariableSalinityBrine::enthalpy( BrineAdapter<FluidState>(fluidState),
                                                                                        VariableSalinityBrine::liquidPhaseIdx )/1e3;

            // mass fraction of CO2 in Brine
            const Scalar X_CO2_w = fluidState.massFraction(liquidPhaseIdx, CO2Idx);

            // heat of dissolution for CO2 according to Fig. 6 in Duan and Sun 2003. (kJ/kg)
            // In the relevant temperature ranges CO2 dissolution is exothermal
            const Scalar delta_hCO2 = (-57.4375 + T * 0.1325) * 1000/44;

            // enthalpy contribution of water and CO2 (kJ/kg)
            const Scalar hw = H2O::liquidEnthalpy(T, p)/1e3;
            const Scalar hg = CO2::liquidEnthalpy(T, p)/1e3 + delta_hCO2;

            // Enthalpy of brine with dissolved CO2 (kJ/kg)
            return (h_ls1 - X_CO2_w*hw + hg*X_CO2_w)*1e3;
        }
        else if (phaseIdx == gasPhaseIdx)
        {
            // we assume NaCl to not enter the gas phase, only consider H2O and CO2
            return H2O::gasEnthalpy(T, p)*fluidState.massFraction(gasPhaseIdx, BrineOrH2OIdx)
                   + CO2::gasEnthalpy(T, p) *fluidState.massFraction(gasPhaseIdx, CO2Idx);
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index.");
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

        if (phaseIdx == liquidPhaseIdx)
        {
            if (componentIdx == BrineOrH2OIdx)
                return H2O::liquidEnthalpy(T, p);
            else if (componentIdx == CO2Idx)
                return CO2::liquidEnthalpy(T, p);
            else if (componentIdx == NaClIdx)
                DUNE_THROW(Dune::NotImplemented, "The component enthalpy for NaCl is not implemented.");
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << componentIdx);
        }
        else if (phaseIdx == gasPhaseIdx)
        {
            if (componentIdx == BrineOrH2OIdx)
                return H2O::gasEnthalpy(T, p);
            else if (componentIdx == CO2Idx)
                return CO2::gasEnthalpy(T, p);
            else if (componentIdx == NaClIdx)
                DUNE_THROW(Dune::InvalidStateException, "Implementation assumes NaCl not to be present in gas phase");
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << componentIdx);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    using Base::thermalConductivity;
    /*!
     * \brief Thermal conductivity of a fluid phase \f$\mathrm{[W/(m K)]}\f$.
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     *
     * \note For the thermal conductivity of the phases the contribution of the minor
     *       component is neglected. This contribution is probably not big, but somebody
     *       would have to find out its influence.
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState& fluidState, int phaseIdx)
    {
        if (phaseIdx == liquidPhaseIdx)
            return useConstantSalinity ? ConstantSalinityBrine::liquidThermalConductivity( fluidState.temperature(phaseIdx),
                                                                                           fluidState.pressure(phaseIdx) )
                                       : VariableSalinityBrine::thermalConductivity( BrineAdapter<FluidState>(fluidState),
                                                                                     VariableSalinityBrine::liquidPhaseIdx );
        else if (phaseIdx == gasPhaseIdx)
            return CO2::gasThermalConductivity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index.");
    }

    using Base::heatCapacity;
    /*!
     * \copybrief Base::heatCapacity
     *
     * \note We employ the heat capacity of the pure phases.
     *
     * \todo TODO Implement heat capacity for gaseous CO2
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        if(phaseIdx == liquidPhaseIdx)
            return useConstantSalinity ? ConstantSalinityBrine::liquidHeatCapacity( fluidState.temperature(phaseIdx),
                                                                                    fluidState.pressure(phaseIdx) )
                                       : VariableSalinityBrine::heatCapacity( BrineAdapter<FluidState>(fluidState),
                                                                              VariableSalinityBrine::liquidPhaseIdx );
        else if (phaseIdx == gasPhaseIdx)
            return CO2::liquidHeatCapacity(fluidState.temperature(phaseIdx),
                                           fluidState.pressure(phaseIdx));

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index.");
    }

private:

    /*!
     * \brief Liquid-phase density calculation of a mixture of brine and CO2 accounting for compositional effects.
     *
     * \param fluidState An arbitrary fluid state
     * \return liquidDensity the liquid-phase density
     */
    template<class FluidState>
    static Scalar liquidDensityMixture_(const FluidState& fluidState)
    {
        const auto T = fluidState.temperature(liquidPhaseIdx);
        const auto p = fluidState.pressure(liquidPhaseIdx);

        if (T < 273.15)
            DUNE_THROW(NumericalProblem, "Liquid density for Brine and CO2 is only "
                                         "defined above 273.15K (T = " << T << ")");

        if (p >= 2.5e8)
            DUNE_THROW(NumericalProblem, "Liquid density for Brine and CO2 is only "
                                         "defined below 250MPa (p = " << p << ")");

        // density of pure water
        Scalar rho_pure = H2O::liquidDensity(T, p);

        // density of water with dissolved CO2 (neglect NaCl)
        Scalar rho_lCO2;
        if (useConstantSalinity)
        {
            // use normalized composition for to calculate the density
            // (the relations don't seem to take non-normalized compositions too well...)
            // TODO Do we really need this normalization???
            using std::min;
            using std::max;
            Scalar xlBrine = min(1.0, max(0.0, fluidState.moleFraction(liquidPhaseIdx, BrineOrH2OIdx)));
            Scalar xlCO2 = min(1.0, max(0.0, fluidState.moleFraction(liquidPhaseIdx, CO2Idx)));
            Scalar sumx = xlBrine + xlCO2;
            xlBrine /= sumx;
            xlCO2 /= sumx;

            rho_lCO2 = liquidDensityWaterCO2_(T, p, xlBrine, xlCO2);
        }
        else
        {
            // rescale mole fractions
            auto xlH2O = fluidState.moleFraction(liquidPhaseIdx, BrineOrH2OIdx);
            auto xlCO2 = fluidState.moleFraction(liquidPhaseIdx, NaClIdx);
            const auto sumMoleFrac = xlH2O + xlCO2;
            xlH2O = xlH2O/sumMoleFrac;
            xlCO2 = xlCO2/sumMoleFrac;
            rho_lCO2 = liquidDensityWaterCO2_(T, p, xlH2O, xlCO2);
        }

        // density of brine (water with nacl)
        Scalar rho_brine = useConstantSalinity ? ConstantSalinityBrine::liquidDensity(T, p)
                                               : VariableSalinityBrine::density( BrineAdapter<FluidState>(fluidState),
                                                                                 VariableSalinityBrine::liquidPhaseIdx );

        // contribution of co2 to the density
        Scalar contribCO2 = rho_lCO2 - rho_pure;

        // Total brine density with dissolved CO2
        // rho_{b, CO2} = rho_pure + contribution(salt) + contribution(CO2)
        return rho_brine + contribCO2;
    }


    /*!
     * \brief Liquid-phase density for a mixture of CO2 in pure water.
     * \note this is used by liquidDensityMixture_
     *
     * \param temperature The temperature
     * \param pl the liquid-phase pressure
     * \param xlH2O the liquid-phase H2O mole fraction
     * \param xlCO2 the liquid-phase CO2 mole fraction
     * \return the density of a mixture of CO2 in pure water
     */
    static Scalar liquidDensityWaterCO2_(Scalar temperature,
                                         Scalar pl,
                                         Scalar xlH2O,
                                         Scalar xlCO2)
    {
        const Scalar M_CO2 = CO2::molarMass();
        const Scalar M_H2O = H2O::molarMass();

        const Scalar tempC = temperature - 273.15;        /* tempC : temperature in °C */
        const Scalar rho_pure = H2O::liquidDensity(temperature, pl);

        // xlH2O is available, but in case of a pure gas phase
        // the value of M_T for the virtual liquid phase can become very large
        xlH2O = 1.0 - xlCO2;
        const Scalar M_T = M_H2O * xlH2O + M_CO2 * xlCO2;
        const Scalar V_phi = (37.51 +
                              tempC*(-9.585e-2 +
                                     tempC*(8.74e-4 -
                                            tempC*5.044e-7))) / 1.0e6;
        return 1/(xlCO2 * V_phi/M_T + M_H2O * xlH2O / (rho_pure * M_T));
    }
};

} // end namespace FluidSystems
} // end namespace Dumux

#endif
