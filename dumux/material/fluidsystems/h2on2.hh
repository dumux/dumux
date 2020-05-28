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
 * \brief @copydoc Dumux::FluidSystems::H2ON2
 */
#ifndef DUMUX_H2O_N2_FLUID_SYSTEM_HH
#define DUMUX_H2O_N2_FLUID_SYSTEM_HH

#include <cassert>
#include <iomanip>

#include <dumux/common/exceptions.hh>

#include <dumux/material/idealgas.hh>

#include <dumux/material/components/n2.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/binarycoefficients/h2o_n2.hh>

#include <dumux/io/name.hh>

#include "base.hh"

namespace Dumux {
namespace FluidSystems {
/*!
 * \ingroup Fluidsystems
 * \brief Policy for the H2O-N2 fluid system
 */
template<bool fastButSimplifiedRelations = false>
struct H2ON2DefaultPolicy
{
    static constexpr  bool useH2ODensityAsLiquidMixtureDensity() { return fastButSimplifiedRelations; }
    static constexpr  bool useIdealGasDensity() { return fastButSimplifiedRelations; }
    static constexpr  bool useN2ViscosityAsGasMixtureViscosity() { return fastButSimplifiedRelations; }
    static constexpr  bool useN2HeatConductivityAsGasMixtureHeatConductivity() { return fastButSimplifiedRelations; }
    static constexpr  bool useIdealGasHeatCapacities() { return fastButSimplifiedRelations; }
};

/*!
 * \ingroup Fluidsystems
 *
 * \brief A two-phase fluid system with two components water \f$(\mathrm{H_2O})\f$
 *        Nitrogen \f$(\mathrm{N_2})\f$ for non-equilibrium models.
 */
template <class Scalar, class Policy = H2ON2DefaultPolicy<>>
class H2ON2
    : public Base<Scalar, H2ON2<Scalar, Policy> >
{
    using ThisType = H2ON2<Scalar, Policy>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;

    // convenience using declarations
    using IdealGas = Dumux::IdealGas<Scalar>;
    using TabulatedH2O = Components::TabulatedComponent<Dumux::Components::H2O<Scalar> >;
    using SimpleN2 = Dumux::Components::N2<Scalar>;

public:
    using H2O = TabulatedH2O; //!< The components for pure water
    using N2 = SimpleN2; //!< The components for pure nitrogen

    static constexpr int numPhases = 2; //!< Number of phases in the fluid system
    static constexpr int numComponents = 2; //!< Number of components in the fluid system

    static constexpr int liquidPhaseIdx = 0; //!< index of the liquid phase
    static constexpr int gasPhaseIdx = 1; //!< index of the gas phase
    static constexpr int phase0Idx = liquidPhaseIdx; //!< index of the first phase
    static constexpr int phase1Idx = gasPhaseIdx; //!< index of the second phase

    static constexpr int H2OIdx = 0;
    static constexpr int N2Idx = 1;
    static constexpr int comp0Idx = H2OIdx; //!< index of the first component
    static constexpr int comp1Idx = N2Idx; //!< index of the second component
    static constexpr int liquidCompIdx = H2OIdx; //!< index of the liquid component
    static constexpr int gasCompIdx = N2Idx; //!< index of the gas component

    /****************************************
     * Fluid phase related static parameters
     ****************************************/
    /*!
     * \brief Return the human readable name of a fluid phase
     *
     * \param phaseIdx The index of the fluid phase to consider
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
    static bool isIdealGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        if (phaseIdx == gasPhaseIdx)
            // let the components decide
            return H2O::gasIsIdeal() && N2::gasIsIdeal();
        return false; // not a gas
    }

    /****************************************
     * Component related static parameters
     ****************************************/
    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        switch (compIdx)
        {
            case H2OIdx: return H2O::name();
            case N2Idx: return N2::name();
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
            H2O::molarMass(),
            N2::molarMass(),
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return M[compIdx];
    }

    /*!
     * \brief Critical temperature of a component \f$\mathrm{[K]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalTemperature(int compIdx)
    {
        static const Scalar Tcrit[] = {
            H2O::criticalTemperature(),
            N2::criticalTemperature()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return Tcrit[compIdx];
    }

    /*!
     * \brief Critical pressure of a component \f$\mathrm{[Pa]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalPressure(int compIdx)
    {
        static const Scalar pcrit[] = {
            H2O::criticalPressure(),
            N2::criticalPressure()
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return pcrit[compIdx];
    }

    /*!
     * \brief Vapor pressure including the Kelvin equation in \f$\mathrm{[Pa]}\f$
     *
     * Calculate the decreased vapor pressure due to capillarity
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar kelvinVaporPressure(const FluidState &fluidState,
                                      const int phaseIdx,
                                      const int compIdx)
    {
        assert(compIdx == H2OIdx && phaseIdx == liquidPhaseIdx);

        using std::exp;
        return fugacityCoefficient(fluidState, phaseIdx, compIdx)
               * fluidState.pressure(phaseIdx)
               * exp(-(fluidState.pressure(gasPhaseIdx)-fluidState.pressure(liquidPhaseIdx))
                       / density(fluidState, phaseIdx)
                       / (Dumux::Constants<Scalar>::R / molarMass(compIdx))
                       / fluidState.temperature());
    }

    /*!
     * \brief Molar volume of a component at the critical point \f$\mathrm{[m^3/mol]}\f$.
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar criticalMolarVolume(int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "H2ON2FluidSystem::criticalMolarVolume()");
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
            N2::acentricFactor()
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
             /*pMin=*/0.0,
             /*pMax=*/20e6,
             /*numP=*/200);
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges
     *
     * \param tempMin The minimum temperature used for tabulation of water \f$\mathrm{[K]}\f$
     * \param tempMax The maximum temperature used for tabulation of water \f$\mathrm{[K]}\f$
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water \f$\mathrm{[Pa]}\f$
     * \param pressMax The maximum pressure used for tabulation of water \f$\mathrm{[Pa]}\f$
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        std::cout << "The H2O-N2 fluid system was configured with the following policy:\n";
        std::cout << " - use H2O density as liquid mixture density: " << std::boolalpha << Policy::useH2ODensityAsLiquidMixtureDensity() << "\n";
        std::cout << " - use ideal gas density: " << std::boolalpha << Policy::useIdealGasDensity() << "\n";
        std::cout << " - use N2 viscosity as gas mixture viscosity: " << std::boolalpha << Policy::useN2ViscosityAsGasMixtureViscosity() << "\n";
        std::cout << " - use N2 heat conductivity as gas mixture heat conductivity: " << std::boolalpha << Policy::useN2HeatConductivityAsGasMixtureHeatConductivity() << "\n";
        std::cout << " - use ideal gas heat capacities: " << std::boolalpha << Policy::useIdealGasHeatCapacities() << std::endl;

        if (H2O::isTabulated)
        {
            TabulatedH2O::init(tempMin, tempMax, nTemp,
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
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == liquidPhaseIdx) {
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
                          + N2::molarMass()*fluidState.moleFraction(liquidPhaseIdx, N2Idx));
            }
        }

        // gas phase
        using std::max;
        if (Policy::useIdealGasDensity())
            // for the gas phase assume an ideal gas
        {
            const Scalar averageMolarMass = fluidState.averageMolarMass(gasPhaseIdx);
            return IdealGas::density(averageMolarMass, T, p);
        }

        // assume ideal mixture: steam and nitrogen don't "see" each other
        Scalar rho_gH2O = H2O::gasDensity(T, fluidState.partialPressure(gasPhaseIdx, H2OIdx));
        Scalar rho_gN2 = N2::gasDensity(T, fluidState.partialPressure(gasPhaseIdx, N2Idx));
        return (rho_gH2O + rho_gN2);
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
        if (phaseIdx == liquidPhaseIdx)
        {
            // assume pure water or that each gas molecule displaces exactly one
            // molecule in the liquid.
            return H2O::liquidMolarDensity(T, p);
        }

        // gas phase
        using std::max;
        if (Policy::useIdealGasDensity())
            // for the gas phase assume an ideal gas
        {
            return IdealGas::molarDensity(T, p);
        }

        // assume ideal mixture: steam and nitrogen don't "see" each other
        Scalar rho_gH2O = H2O::gasMolarDensity(T, fluidState.partialPressure(gasPhaseIdx, H2OIdx));
        Scalar rho_gN2 = N2::gasMolarDensity(T, fluidState.partialPressure(gasPhaseIdx, N2Idx));
        return rho_gH2O + rho_gN2;
    }

    using Base::viscosity;
    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase \f$\mathrm{[Pa*s]}\f$
     *
     * Compositional effects in the gas phase are accounted by the Wilke method.
     * See Reid et al. (1987)  \cite reid1987 <BR>
     * 4th edition, McGraw-Hill, 1987, 407-410
     * 5th edition, McGraw-Hill, 20001, p. 9.21/22
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

        // liquid phase
        if (phaseIdx == liquidPhaseIdx) {
            // assume pure water for the liquid phase
            return H2O::liquidViscosity(T, p);
        }

        // gas phase
        if (Policy::useN2ViscosityAsGasMixtureViscosity())
        {
            // assume pure nitrogen for the gas phase
            return N2::gasViscosity(T, p);
        }
        else
        {
            // Wilke method (Reid et al.):
            Scalar muResult = 0;
            const Scalar mu[numComponents] = {
                H2O::gasViscosity(T, H2O::vaporPressure(T)),
                N2::gasViscosity(T, p)
            };

            Scalar sumx = 0.0;
            using std::max;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                sumx += fluidState.moleFraction(phaseIdx, compIdx);
            sumx = max(1e-10, sumx);

            for (int i = 0; i < numComponents; ++i) {
                Scalar divisor = 0;
//		using std::sqrt;
//		using std::pow;
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
    }

    using Base::fugacityCoefficient;
    /*!
     * \brief Calculate the fugacity coefficient \f$\mathrm{[-]}\f$ of an individual
     *        component in a fluid phase
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
     * The quantity "fugacity" itself is just an other way to express
     * the chemical potential \f$\mathrm{\zeta^\kappa_\alpha}\f$ of the
     * component. It is defined via
     *
     * \f[
     f^\kappa_\alpha := \exp\left\{\frac{\zeta^\kappa_\alpha}{k_B T_\alpha} \right\}
     \f]
     * where \f$\mathrm{k_B = 1.380\cdot10^{-23}\;J/K}\f$ is the Boltzmann constant.
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
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == liquidPhaseIdx) {
            if (compIdx == H2OIdx)
                return H2O::vaporPressure(T)/p;
            return BinaryCoeff::H2O_N2::henry(T)/p;
        }

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        return 1.0;
    }

    using Base::diffusionCoefficient;
    /*!
     * \brief Calculate the molecular diffusion coefficient for a
     *        component in a fluid phase \f$\mathrm{[mol^2 * s / (kg*m^3)]}\f$
     *
     * Molecular diffusion of a compoent \f$\mathrm{\kappa}\f$ is caused by a
     * gradient of the chemical potential and follows the law
     *
     * \f[ J = - D \mathbf{grad} \mu_\kappa \f]
     *
     * where \f$\mathrm{\mu_\kappa}\f$ is the component's chemical potential,
     * \f$\mathrm{D}\f$ is the diffusion coefficient and \f$\mathrm{J}\f$ is the
     * diffusive flux. \f$\mathrm{mu_\kappa}\f$ is connected to the component's
     * fugacity \f$\mathrm{f_\kappa}\f$ by the relation
     *
     * \f[ \mu_\kappa = R T_\alpha \mathrm{ln} \frac{f_\kappa}{p_\alpha} \f]
     *
     * where \f$p_\alpha\f$ and \f$T_\alpha\f$ are the fluid phase'
     * pressure and temperature.
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
        DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients");
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
        if (compIIdx > compJIdx)
        {
            using std::swap;
            swap(compIIdx, compJIdx);
        }

        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == liquidPhaseIdx && compIIdx == H2OIdx && compJIdx == N2Idx)
            return BinaryCoeff::H2O_N2::liquidDiffCoeff(T, p);

        else if (phaseIdx == gasPhaseIdx && compIIdx == H2OIdx && compJIdx == N2Idx)
            return BinaryCoeff::H2O_N2::gasDiffCoeff(T, p);

        else
            DUNE_THROW(Dune::InvalidStateException,
                   "Binary diffusion coefficient of components "
                    << compIIdx << " and " << compJIdx
                    << " in phase " << phaseIdx << " is unavailable!\n");
    }

    using Base::enthalpy;
    /*!
     * \brief Given a phase's composition, temperature, pressure and
     *        density, calculate its specific enthalpy \f$\mathrm{[J/kg]}\f$.
     *
     *  \note This fluid system neglects the contribution of
     *        gas-molecules in the liquid phase. This contribution is
     *        probably not big. Somebody would have to find out the
     *        enthalpy of solution for this system. ...
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           int phaseIdx)
    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == liquidPhaseIdx) {
            return H2O::liquidEnthalpy(T, p);
        }
        // gas phase
        else {
            // assume ideal mixture: which means
            // that the total specific enthalpy is the sum of the
            // "partial specific enthalpies" of the components.
            Scalar hH2O =
                fluidState.massFraction(gasPhaseIdx, H2OIdx)
                * H2O::gasEnthalpy(T, p);
            Scalar hN2 =
                fluidState.massFraction(gasPhaseIdx, N2Idx)
                * N2::gasEnthalpy(T, p);
            return hH2O + hN2;
        }
    }

    /*!
    * \brief Returns the specific enthalpy \f$\mathrm{[J/kg]}\f$ of a component in the specified phase
    * \param fluidState The fluid state
    * \param phaseIdx The index of the phase
    * \param componentIdx The index of the component
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
            if (componentIdx == H2OIdx)
                return H2O::liquidEnthalpy(T, p);
            else if (componentIdx == N2Idx)
                DUNE_THROW(Dune::NotImplemented, "Component enthalpy of nitrogen in liquid phase");
            else
                DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << componentIdx);
        }
        else if (phaseIdx == gasPhaseIdx)
        {
            if (componentIdx == H2OIdx)
                return H2O::gasEnthalpy(T, p);
            else if (componentIdx == N2Idx)
                return N2::gasEnthalpy(T, p);
            else
                DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << componentIdx);
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    using Base::thermalConductivity;
    /*!
     * \brief Thermal conductivity of a fluid phase \f$\mathrm{[W/(m K)]}\f$.
     *
     * Use the conductivity of air and water as a first approximation.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        const Scalar temperature  = fluidState.temperature(phaseIdx) ;
        const Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == liquidPhaseIdx)
        {
            return H2O::liquidThermalConductivity(temperature, pressure);
        }
        else
        {
            Scalar lambdaPureN2 = N2::gasThermalConductivity(temperature, pressure);
            if (!Policy::useN2HeatConductivityAsGasMixtureHeatConductivity())
            {
                Scalar xN2 = fluidState.moleFraction(phaseIdx, N2Idx);
                Scalar xH2O = fluidState.moleFraction(phaseIdx, H2OIdx);
                Scalar lambdaN2 = xN2 * lambdaPureN2;
                Scalar partialPressure  = pressure * xH2O;
                Scalar lambdaH2O = xH2O * H2O::gasThermalConductivity(temperature, partialPressure);
                return lambdaN2 + lambdaH2O;
            }
            else
                return lambdaPureN2;
        }
    }

    using Base::heatCapacity;
    /*!
     * \brief Specific isobaric heat capacity of a fluid phase.
     *        \f$\mathrm{[J/(kg K)]}\f$.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        if (phaseIdx == liquidPhaseIdx) {
            return H2O::liquidHeatCapacity(fluidState.temperature(phaseIdx),
                                           fluidState.pressure(phaseIdx));
        }

        // for the gas phase, assume ideal mixture
        Scalar c_pN2;
        Scalar c_pH2O;
        // let the water and nitrogen components do things their own way
        if (!Policy::useIdealGasHeatCapacities()) {
            c_pN2 = N2::gasHeatCapacity(fluidState.temperature(phaseIdx),
                                        fluidState.pressure(phaseIdx)
                                        * fluidState.moleFraction(phaseIdx, N2Idx));

            c_pH2O = H2O::gasHeatCapacity(fluidState.temperature(phaseIdx),
                                          fluidState.pressure(phaseIdx)
                                          * fluidState.moleFraction(phaseIdx, H2OIdx));
        }
        else {
            // assume an ideal gas for both components. See:
            // http://en.wikipedia.org/wiki/Heat_capacity
            Scalar c_vN2molar = Constants<Scalar>::R*2.39;
            Scalar c_pN2molar = Constants<Scalar>::R + c_vN2molar;

            Scalar c_vH2Omolar = Constants<Scalar>::R*3.37; // <- correct??
            Scalar c_pH2Omolar = Constants<Scalar>::R + c_vH2Omolar;

            c_pN2 = c_pN2molar/molarMass(N2Idx);
            c_pH2O = c_pH2Omolar/molarMass(H2OIdx);
        }

        // mangle both components together
        return c_pH2O*fluidState.massFraction(gasPhaseIdx, H2OIdx)
               + c_pN2*fluidState.massFraction(gasPhaseIdx, N2Idx);
    }
};

} // end namespace FluidSystems

} // end namespace Dumux

#endif
