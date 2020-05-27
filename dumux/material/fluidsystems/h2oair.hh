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
 * \brief @copybrief Dumux::FluidSystems::H2OAir
 */
#ifndef DUMUX_H2O_AIR_SYSTEM_HH
#define DUMUX_H2O_AIR_SYSTEM_HH

#include <cassert>
#include <iomanip>

#include <dumux/material/idealgas.hh>
#include <dumux/material/fluidsystems/base.hh>

#include <dumux/material/binarycoefficients/h2o_air.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/h2o.hh>

#include <dumux/common/exceptions.hh>

#include <dumux/io/name.hh>

namespace Dumux {
namespace FluidSystems {
/*!
 * \ingroup Fluidsytems
 * \brief Policy for the H2O-air fluid system
 */
template<bool fastButSimplifiedRelations = false>
struct H2OAirDefaultPolicy
{
    static constexpr bool useH2ODensityAsLiquidMixtureDensity() { return fastButSimplifiedRelations; }
    static constexpr bool useIdealGasDensity() { return fastButSimplifiedRelations; }
    static constexpr bool useAirViscosityAsGasMixtureViscosity() { return fastButSimplifiedRelations; }
};

/*!
 * \ingroup Fluidsystems
 *
 * \brief A compositional two-phase fluid system with water and air as
 *        components in both, the liquid and the gas phase.
 *
 * This fluidsystem features gas and liquid phases of distilled water
 * \f$(\mathrm{H_2O})\f$) and air (Pseudo component composed of \f$\mathrm{79\%\;N_2}\f$,
 * \f$\mathrm{20\%\;O_2}\f$ and \f$\mathrm{1\%\;Ar}\f$) as components. It is applied by
 * default with the tabulated version of water of the IAPWS-formulation.
 */
template <class Scalar,
          class H2Otype = Components::TabulatedComponent<Components::H2O<Scalar> >,
          class Policy = H2OAirDefaultPolicy<>,
          bool useKelvinVaporPressure = false>
class H2OAir
: public Base<Scalar, H2OAir<Scalar, H2Otype, Policy> >
{
    using ThisType = H2OAir<Scalar,H2Otype, Policy>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;
    using IdealGas = Dumux::IdealGas<Scalar>;

public:
    using H2O = H2Otype;
    using Air = Dumux::Components::Air<Scalar>;

    static constexpr int numPhases = 2; //!< Number of phases in the fluid system
    static constexpr int numComponents = 2; //!< Number of components in the fluid system

    static constexpr int liquidPhaseIdx = 0; //!< index of the liquid phase
    static constexpr int gasPhaseIdx = 1; //!< index of the gas phase
    static constexpr int phase0Idx = liquidPhaseIdx; //!< index of the first phase
    static constexpr int phase1Idx = gasPhaseIdx; //!< index of the second phase

    static constexpr int H2OIdx = 0; //!< index of the frist component
    static constexpr int AirIdx = 1; //!< index of the second component
    static constexpr int comp0Idx = H2OIdx; //!< index of the frist component
    static constexpr int comp1Idx = AirIdx; //!< index of the second component
    static constexpr int liquidCompIdx = H2OIdx; //!< index of the liquid component
    static constexpr int gasCompIdx = AirIdx; //!< index of the gas component

    /*!
     * \brief Return the human readable name of a phase
     *
     * \param phaseIdx index of the phase
     */
    static std::string phaseName(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
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
    static constexpr bool isIdealMixture(int phaseIdx)
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
        // ideal gases are always compressible
        if (phaseIdx == gasPhaseIdx)
            return true;
        // the water component decides for the liquid phase...
        return H2O::liquidIsCompressible();
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal gas.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isIdealGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        // let the fluids decide
        if (phaseIdx == gasPhaseIdx)
            return H2O::gasIsIdeal() && Air::gasIsIdeal();
        return false; // not a gas
    }

    /****************************************
     * Component related static parameters
     ****************************************/
    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx index of the component
     */
    static std::string componentName(int compIdx)
    {
        switch (compIdx)
        {
            case H2OIdx: return H2O::name();
            case AirIdx: return Air::name();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the molar mass of a component \f$\mathrm{[kg/mol]}\f$.
     *
     * \param compIdx index of the component
     */
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx)
        {
            case H2OIdx: return H2O::molarMass();
            case AirIdx: return Air::molarMass();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Critical temperature of a component \f$\mathrm{[K]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalTemperature(int compIdx)
    {
        static const Scalar TCrit[] = {
            H2O::criticalTemperature(),
            Air::criticalTemperature()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return TCrit[compIdx];
    }

    /*!
     * \brief Critical pressure of a component \f$\mathrm{[Pa]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalPressure(int compIdx)
    {
        static const Scalar pCrit[] = {
            H2O::criticalPressure(),
            Air::criticalPressure()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return pCrit[compIdx];
    }

    /*!
     * \brief Vapor pressure of a component \f$\mathrm{[Pa]}\f$.
     *
     * \param fluidState The fluid state
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar vaporPressure(const FluidState &fluidState, int compIdx)
    {
        if (compIdx == H2OIdx)
        {
            const auto t = fluidState.temperature(H2OIdx);
            if (!useKelvinVaporPressure)
                return H2O::vaporPressure(t);
            else
            {
                const auto pc =  (fluidState.wettingPhase() == (int) H2OIdx)
                                 ? fluidState.pressure(AirIdx)-fluidState.pressure(H2OIdx)
                                 : fluidState.pressure(H2OIdx)-fluidState.pressure(AirIdx);
                return H2O::vaporPressure(t)*exp( -pc * molarMass(H2OIdx)
                                                      / density(fluidState, H2OIdx)
                                                      / (Dumux::Constants<Scalar>::R*t) );
            }
        }
        else if (compIdx == AirIdx)
            // return Air::vaporPressure(fluidState.temperature(AirIdx));
            DUNE_THROW(Dune::NotImplemented, "Air::vaporPressure(t)");
        else
            DUNE_THROW(Dune::NotImplemented, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Molar volume of a component at the critical point \f$\mathrm{[m^3/mol]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalMolarVolume(int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "H2OAirFluidSystem::criticalMolarVolume()");
    }

    /*!
     * \brief The acentric factor of a component \f$\mathrm{[-]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar acentricFactor(int compIdx)
    {
        static const Scalar accFac[] = {
            H2O::acentricFactor(),
            Air::acentricFactor()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return accFac[compIdx];
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

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
             /*pMin=*/-10.,
             /*pMax=*/20e6,
             /*numP=*/200);
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges
     *
     * \param tempMin The minimum temperature used for tabulation of water \f$\mathrm{[K]}\f$
     * \param tempMax The maximum temperature used for tabulation of water\f$\mathrm{[K]}\f$
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water \f$\mathrm{[Pa]}\f$
     * \param pressMax The maximum pressure used for tabulation of water \f$\mathrm{[Pa]}\f$
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        std::cout << "The H2O-air fluid system was configured with the following policy:\n";
        std::cout << " - use H2O density as liquid mixture density: " << std::boolalpha << Policy::useH2ODensityAsLiquidMixtureDensity() << "\n";
        std::cout << " - use ideal gas density: " << std::boolalpha << Policy::useIdealGasDensity() << "\n";
        std::cout << " - use air viscosity as gas mixture viscosity: " << std::boolalpha << Policy::useAirViscosityAsGasMixtureViscosity() << std::endl;

        if (H2O::isTabulated)
        {
            H2O::init(tempMin, tempMax, nTemp,
                      pressMin, pressMax, nPress);
        }
    }

    using Base::density;
    /*!
     * \brief Given a phase's composition, temperature, pressure, and
     *        the partial pressures of all components, return its
     *        density \f$\mathrm{[kg/m^3]}\f$.
     *
     * If Policy::useH2ODensityAsLiquidMixtureDensity() == false, we apply Eq. (7)
     * in Class et al. (2002a) \cite A3:class:2002b <BR>
     * for the liquid density.
     *
     * \param phaseIdx index of the phase
     * \param fluidState the fluid state
     *
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == phase0Idx)
        {
            if (Policy::useH2ODensityAsLiquidMixtureDensity())
                // assume pure water
                return H2O::liquidDensity(T, p);
            else
            {
                // See: Eq. (7) in Class et al. (2002a)
                // This assumes each gas molecule displaces exactly one
                // molecule in the liquid.
                return H2O::liquidMolarDensity(T, p)
                       * (H2O::molarMass()*fluidState.moleFraction(liquidPhaseIdx, H2OIdx)
                          + Air::molarMass()*fluidState.moleFraction(liquidPhaseIdx, AirIdx));
            }
        }
        else if (phaseIdx == gasPhaseIdx)
        {
            if (Policy::useIdealGasDensity())
                // for the gas phase assume an ideal gas
            {
                const Scalar averageMolarMass = fluidState.averageMolarMass(gasPhaseIdx);
                return IdealGas::density(averageMolarMass, T, p);
            }

            return H2O::gasDensity(T, fluidState.partialPressure(gasPhaseIdx, H2OIdx))
                   + Air::gasDensity(T, fluidState.partialPressure(gasPhaseIdx, AirIdx));
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    using Base::molarDensity;
    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density for the simple relation is defined by the
     * mass density \f$\rho_\alpha\f$ and the molar mass of the main component \f$M_\kappa\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{M_\kappa} \;.\f]
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == phase0Idx)
        {
            // assume pure water or that each gas molecule displaces exactly one
            // molecule in the liquid.
            return H2O::liquidMolarDensity(T, p);
        }
        else if (phaseIdx == phase1Idx)
        {
            if (Policy::useIdealGasDensity())
                // for the gas phase assume an ideal gas
            { return IdealGas::molarDensity(T, p); }

            return H2O::gasMolarDensity(T, fluidState.partialPressure(phase1Idx, H2OIdx))
                   + Air::gasMolarDensity(T, fluidState.partialPressure(phase1Idx, AirIdx));
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    using Base::viscosity;
    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase \f$\mathrm{[Pa*s]}\f$
     *
     * Compositional effects in the gas phase are accounted by the Wilke method.
     * See Reid et al. (1987)  \cite reid1987 <BR>
     * 4th edition, McGraw-Hill, 1987, 407-410
     * 5th edition, McGraw-Hill, 2001, p. 9.21/22
     * \note Compositional effects for a liquid mixture have to be implemented.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == liquidPhaseIdx)
        {
            // assume pure water for the liquid phase
            return H2O::liquidViscosity(T, p);
        }
        else if (phaseIdx == gasPhaseIdx)
        {
            if(Policy::useAirViscosityAsGasMixtureViscosity()){
                return Air::gasViscosity(T, p);
            }
            else //using a complicated version of this fluid system
            {
                // Wilke method (Reid et al.):
                Scalar muResult = 0;
                const Scalar mu[numComponents] = {
                    H2O::gasViscosity(T,
                                      H2O::vaporPressure(T)),
                    Air::gasViscosity(T, p)
                };

                // molar masses
                const Scalar M[numComponents] =  {
                    H2O::molarMass(),
                    Air::molarMass()
                };

                for (int i = 0; i < numComponents; ++i)
                {
                    Scalar divisor = 0;
                    using std::sqrt;
                    for (int j = 0; j < numComponents; ++j)
                    {
                        // 1 + (mu[i]/mu[j]^1/2 * (M[i]/M[j])^1/4)
                        Scalar phiIJ = 1 + sqrt(mu[i]/mu[j] * sqrt(M[j]/M[i]));
                        phiIJ = phiIJ * phiIJ / sqrt(8*(1 + M[i]/M[j]));
                        divisor += fluidState.moleFraction(phaseIdx, j)*phiIJ;
                    }
                    muResult += fluidState.moleFraction(phaseIdx, i)*mu[i] / divisor;
                }
                return muResult;
            }
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    using Base::fugacityCoefficient;
    /*!
     * \brief Returns the fugacity coefficient \f$\mathrm{[-]}\f$ of a component in a
     *        phase.
     *
     * The fugacity coefficient \f$\phi^\kappa_\alpha\f$ of
     * component \f$\kappa\f$ in phase \f$\alpha\f$ is connected to
     * the fugacity \f$f^\kappa_\alpha\f$ and the component's mole
     * fraction \f$x^\kappa_\alpha\f$ by means of the relation
     *
     * \f[
     f^\kappa_\alpha = \phi^\kappa_\alpha\;x^\kappa_\alpha\;p_\alpha
     \f]
     * where \f$p_\alpha\f$ is the pressure of the fluid phase.
     *
     * For liquids with very low miscibility this boils down to the
     * Henry constant for the solutes and the saturated vapor pressure
     * both divided by phase pressure.
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

        if (phaseIdx == liquidPhaseIdx) {
            if (compIdx == H2OIdx)
                return vaporPressure(fluidState, compIdx)/p;
            return BinaryCoeff::H2O_Air::henry(T)/p;
        }

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        return 1.0;
    }

    /*!
     * \brief Returns the relative humidity of the gas phase.
     *
     * The relative humidity is the ratio of the partial pressure of water vapor
     * to the equilibrium vapor pressure of water at a given temperature.
     */
    template <class FluidState>
    static Scalar relativeHumidity(const FluidState &fluidState)
    {
        return fluidState.partialPressure(gasPhaseIdx, comp0Idx)
               / H2O::vaporPressure(fluidState.temperature(gasPhaseIdx));
    }

    using Base::diffusionCoefficient;
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::H2OAir::diffusionCoefficient()");
    }

    using Base::binaryDiffusionCoefficient;
    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for components
     *        \f$i\f$ and \f$j\f$ in this phase.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIIdx The index of the first component to consider
     * \param compJIdx The index of the second component to consider
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        using std::swap;
        if (compIIdx > compJIdx)
            swap(compIIdx, compJIdx);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        // we are in the liquid phase
        if (phaseIdx == liquidPhaseIdx)
        {
            if (compIIdx == H2OIdx && compJIdx == AirIdx)
                return BinaryCoeff::H2O_Air::liquidDiffCoeff(T, p);
            else
                DUNE_THROW(Dune::InvalidStateException,
                           "Binary diffusion coefficient of components "
                            << compIIdx << " and " << compJIdx
                            << " in phase " << phaseIdx << " is undefined!\n");
        }

        // we are in the gas phase
        else if (phaseIdx == gasPhaseIdx)
        {
            if (compIIdx == H2OIdx && compJIdx == AirIdx)
                return BinaryCoeff::H2O_Air::gasDiffCoeff(T, p);
            else
                DUNE_THROW(Dune::InvalidStateException,
                           "Binary diffusion coefficient of components "
                           << compIIdx << " and " << compJIdx
                           << " in phase " << phaseIdx << " is undefined!\n");
        }

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    using Base::enthalpy;
    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return its specific enthalpy \f$\mathrm{[J/kg]}\f$.
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     *
     * See:
     * Class 2001:
     * Theorie und numerische Modellierung nichtisothermer Mehrphasenprozesse in NAPL-kontaminierten porösen Medien
     * Chapter 2.1.13 Innere Energie, Wäremekapazität, Enthalpie \cite A3:class:2001 <BR>
     *
     * Formula (2.42):
     * the specific enthalpy of a gasphase result from the sum of (enthalpies*mass fraction) of the components
     *
     *  \todo This system neglects the contribution of gas-molecules in the liquid phase.
     *        This contribution is probably not big. Somebody would have to find out the enthalpy of solution for this system. ...
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           int phaseIdx)
    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == liquidPhaseIdx)
            return H2O::liquidEnthalpy(T, p);

        else if (phaseIdx == gasPhaseIdx)
            return H2O::gasEnthalpy(T, p)*fluidState.massFraction(gasPhaseIdx, H2OIdx)
                   + Air::gasEnthalpy(T, p)*fluidState.massFraction(gasPhaseIdx, AirIdx);

        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
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

        if (phaseIdx == liquidPhaseIdx)
        {
            // the liquid enthalpy is constant
            return H2O::liquidEnthalpy(T, p);
        }
        else if (phaseIdx == gasPhaseIdx)
        {
            if (componentIdx == H2OIdx)
            {
                return H2O::gasEnthalpy(T, p);
            }
            else if (componentIdx == AirIdx)
            {
                return Air::gasEnthalpy(T, p);
            }
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
     * Use the conductivity of air and water as a first approximation.
     * Source:
     * http://en.wikipedia.org/wiki/List_of_thermal_conductivities
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        const Scalar temperature  = fluidState.temperature(phaseIdx) ;
        const Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == liquidPhaseIdx)
        {
            return H2O::liquidThermalConductivity(temperature, pressure);
        }
        else if (phaseIdx == gasPhaseIdx)
        {
            return Air::gasThermalConductivity(temperature, pressure);
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    using Base::heatCapacity;
    /*!
     * \brief Specific isobaric heat capacity of a fluid phase.
     *        \f$\mathrm{[J/(kg*K)}\f$.
     *
     * \todo Check whether the gas phase enthalpy is a linear mixture of the component
     *       enthalpies and the mole fractions is a good assumption.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx  for which phase to give back the heat capacity
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        const Scalar temperature  = fluidState.temperature(phaseIdx);
        const Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == liquidPhaseIdx)
        {
            // influence of air is neglected
            return H2O::liquidHeatCapacity(temperature, pressure);
        }
        else if (phaseIdx == gasPhaseIdx)
        {
            return Air::gasHeatCapacity(temperature, pressure) * fluidState.moleFraction(gasPhaseIdx, AirIdx)
                   + H2O::gasHeatCapacity(temperature, pressure) * fluidState.moleFraction(gasPhaseIdx, H2OIdx);
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }
};

} // end namespace FluidSystems
} // end namespace Dumux

#endif
