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
 * \brief @copybrief Dumux::FluidSystems::H2OAirXylene
 */
#ifndef DUMUX_H2O_AIR_XYLENE_FLUID_SYSTEM_HH
#define DUMUX_H2O_AIR_XYLENE_FLUID_SYSTEM_HH

#include <dumux/material/idealgas.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/xylene.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

#include <dumux/material/binarycoefficients/h2o_air.hh>
#include <dumux/material/binarycoefficients/h2o_xylene.hh>
#include <dumux/material/binarycoefficients/air_xylene.hh>

#include <dumux/material/fluidsystems/base.hh>

#include <dumux/io/name.hh>

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 * \brief A three-phase fluid system featuring gas, NAPL and water as phases and
 *        distilled water \f$(\mathrm{H_2O})\f$ and air (Pseudo component composed of
 *        \f$\mathrm{79\%\;N_2}\f$, \f$\mathrm{20\%\;O_2}\f$ and Mesitylene \f$(\mathrm{C_8H_{10}})\f$ as components.
 *
 * \note This fluid system assumes all phases to be ideal mixtures.
 */
template <class Scalar,
          class H2OType = Components::TabulatedComponent<Components::H2O<Scalar> > >
class H2OAirXylene
    : public Base<Scalar, H2OAirXylene<Scalar, H2OType> >
{
    using ThisType = H2OAirXylene<Scalar, H2OType>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;

public:
    using H2O = H2OType;
    using NAPL = Components::Xylene<Scalar>;
    using Air = Dumux::Components::Air<Scalar>;

    static const int numPhases = 3;
    static const int numComponents = 3;

    static const int wPhaseIdx = 0; // index of the water phase
    static const int nPhaseIdx = 1; // index of the NAPL phase
    static const int gPhaseIdx = 2; // index of the gas phase

    static const int H2OIdx = 0;
    static const int NAPLIdx = 1;
    static const int AirIdx = 2;

    // export component indices to indicate the main component
    // of the corresponding phase at atmospheric pressure 1 bar
    // and room temperature 20°C:
    static const int wCompIdx = H2OIdx;
    static const int nCompIdx = NAPLIdx;
    static const int gCompIdx = AirIdx;

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
        if (H2O::isTabulated)
        {
            H2O::init(tempMin, tempMax, nTemp,
                      pressMin, pressMax, nPress);
        }
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
    { return phaseIdx == gPhaseIdx && H2O::gasIsIdeal() && Air::gasIsIdeal() && NAPL::gasIsIdeal(); }

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
        if (phaseIdx == gPhaseIdx)
            return true;
        else if (phaseIdx == wPhaseIdx)
            // the water component decides for the water phase...
            return H2O::liquidIsCompressible();

        // the NAPL component decides for the napl phase...
        return NAPL::liquidIsCompressible();
    }

    /*!
     * \brief Return the human readable name of a phase (used in indices)
     * \param phaseIdx The index of the fluid phase to consider
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
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        switch (compIdx) {
            case H2OIdx: return H2O::name();
            case AirIdx: return Air::name();
            case NAPLIdx: return NAPL::name();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the molar mass of a component in \f$\mathrm{[kg/mol]}\f$
     * \param compIdx The index of the component to consider
     */
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx) {
            case H2OIdx: return H2O::molarMass();
            case AirIdx: return Air::molarMass();
            case NAPLIdx: return NAPL::molarMass();
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    using Base::density;
    /*!
     * \brief Given a phase's composition, temperature, pressure, and
     *        the partial pressures of all components, return its
     *        density \f$\mathrm{[kg/m^3]}\f$.
     *
     * We apply Eq. (7)
     * in Class et al. (2002a) \cite A3:class:2002b <BR>
     * for the water density.
     *
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState, int phaseIdx)
    {
        if (phaseIdx == wPhaseIdx) {
            // This assumes each gas molecule displaces exactly one
            // molecule in the liquid.
            return H2O::liquidMolarDensity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx))
                   * (H2O::molarMass()*fluidState.moleFraction(wPhaseIdx, H2OIdx)
                     + Air::molarMass()*fluidState.moleFraction(wPhaseIdx, AirIdx)
                     + NAPL::molarMass()*fluidState.moleFraction(wPhaseIdx, NAPLIdx));
        }
        else if (phaseIdx == nPhaseIdx) {
            // assume pure NAPL for the NAPL phase
            Scalar pressure = NAPL::liquidIsCompressible()?fluidState.pressure(phaseIdx):1e100;
            return NAPL::liquidDensity(fluidState.temperature(phaseIdx), pressure);
        }

        assert (phaseIdx == gPhaseIdx);
        Scalar pH2O =
            fluidState.moleFraction(gPhaseIdx, H2OIdx)  *
            fluidState.pressure(gPhaseIdx);
        Scalar pAir =
            fluidState.moleFraction(gPhaseIdx, AirIdx)  *
            fluidState.pressure(gPhaseIdx);
        Scalar pNAPL =
            fluidState.moleFraction(gPhaseIdx, NAPLIdx)  *
            fluidState.pressure(gPhaseIdx);
        return H2O::gasDensity(fluidState.temperature(phaseIdx), pH2O)
               + Air::gasDensity(fluidState.temperature(phaseIdx), pAir)
               + NAPL::gasDensity(fluidState.temperature(phaseIdx), pNAPL);
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
            return NAPL::liquidMolarDensity(temperature, pressure);
        }
        else if (phaseIdx == wPhaseIdx)
        {   // This assumes each gas molecule displaces exactly one
            // molecule in the liquid.
            return H2O::liquidMolarDensity(temperature, pressure);
        }
        else
        {
            return H2O::gasMolarDensity(temperature, fluidState.partialPressure(gPhaseIdx, H2OIdx))
                   + NAPL::gasMolarDensity(temperature, fluidState.partialPressure(gPhaseIdx, NAPLIdx))
                   + Air::gasMolarDensity(temperature, fluidState.partialPressure(gPhaseIdx, AirIdx));
        }
    }

    using Base::viscosity;
    /*!
     * \brief Return the viscosity of a phase \f$\mathrm{[Pa s]}\f$.
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase to consider
     *
     * \todo Check the parameter phiCAW for the xylene case and give a physical meaningful name
     */
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
            return NAPL::liquidViscosity(fluidState.temperature(phaseIdx),
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
        Scalar muResult;
        const Scalar mu[numComponents] = {
            H2O::gasViscosity(fluidState.temperature(phaseIdx), H2O::vaporPressure(fluidState.temperature(phaseIdx))),
            Air::gasViscosity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx)),
            NAPL::gasViscosity(fluidState.temperature(phaseIdx), NAPL::vaporPressure(fluidState.temperature(phaseIdx)))
        };
        // molar masses
        const Scalar M[numComponents] = {
            H2O::molarMass(),
            Air::molarMass(),
            NAPL::molarMass()
        };

        Scalar muAW = mu[AirIdx]*fluidState.moleFraction(gPhaseIdx, AirIdx)
                      + mu[H2OIdx]*fluidState.moleFraction(gPhaseIdx, H2OIdx)
                        / (fluidState.moleFraction(gPhaseIdx, AirIdx)
                      + fluidState.moleFraction(gPhaseIdx, H2OIdx));
        Scalar xAW = fluidState.moleFraction(gPhaseIdx, AirIdx)
                     + fluidState.moleFraction(gPhaseIdx, H2OIdx);

        Scalar MAW = (fluidState.moleFraction(gPhaseIdx, AirIdx)*Air::molarMass()
                      + fluidState.moleFraction(gPhaseIdx, H2OIdx)*H2O::molarMass())
                     / xAW;

        Scalar phiCAW = 0.3; // simplification for this particular system
        /* actually like this
        * using std::sqrt;
        * using std::pow;
         * Scalar phiCAW = pow(1.+sqrt(mu[NAPLIdx]/muAW)*pow(MAW/M[NAPLIdx],0.25),2)
         *                 / sqrt(8.*(1.+M[NAPLIdx]/MAW));
         */
        Scalar phiAWC = phiCAW * muAW*M[NAPLIdx]/(mu[NAPLIdx]*MAW);

        muResult = (xAW*muAW)/(xAW+fluidState.moleFraction(gPhaseIdx, NAPLIdx)*phiAWC)
                   + (fluidState.moleFraction(gPhaseIdx, NAPLIdx) * mu[NAPLIdx])
                     / (fluidState.moleFraction(gPhaseIdx, NAPLIdx) + xAW*phiCAW);
        return muResult;
    }


    using Base::diffusionCoefficient;
    /*!
     * \brief Given all mole fractions, return the diffusion
     *        coefficient \f$\mathrm{[m^2/s]}\f$ of a component in a phase.
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        Scalar diffCont;

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx==gPhaseIdx) {
            Scalar diffAC = BinaryCoeff::Air_Xylene::gasDiffCoeff(temperature, pressure);
            Scalar diffWC = BinaryCoeff::H2O_Xylene::gasDiffCoeff(temperature, pressure);
            Scalar diffAW = BinaryCoeff::H2O_Air::gasDiffCoeff(temperature, pressure);

            const Scalar xga = fluidState.moleFraction(gPhaseIdx, AirIdx);
            const Scalar xgw = fluidState.moleFraction(gPhaseIdx, H2OIdx);
            const Scalar xgc = fluidState.moleFraction(gPhaseIdx, NAPLIdx);

            if (compIdx==NAPLIdx)
                return (1.- xgw)/(xga/diffAW + xgc/diffWC);
            else if (compIdx==H2OIdx)
                return (1.- xgc)/(xgw/diffWC + xga/diffAC);
            else if (compIdx==AirIdx)
                DUNE_THROW(Dune::InvalidStateException,
                                                 "Diffusivity of Air in the gas phase "
                                                 "is constraint by sum of diffusive fluxes = 0 !\n");
        } else if (phaseIdx==wPhaseIdx){
            Scalar diffACl = BinaryCoeff::Air_Xylene::liquidDiffCoeff(temperature, pressure);
            Scalar diffWCl = BinaryCoeff::H2O_Xylene::liquidDiffCoeff(temperature, pressure);
            Scalar diffAWl = BinaryCoeff::H2O_Air::liquidDiffCoeff(temperature, pressure);

            Scalar xwa = fluidState.moleFraction(wPhaseIdx, AirIdx);
            Scalar xww = fluidState.moleFraction(wPhaseIdx, H2OIdx);
            Scalar xwc = fluidState.moleFraction(wPhaseIdx, NAPLIdx);

            switch (compIdx) {
            case NAPLIdx:
                diffCont = (1.- xww)/(xwa/diffAWl + xwc/diffWCl);
                return diffCont;
            case AirIdx:
                diffCont = (1.- xwc)/(xww/diffWCl + xwa/diffACl);
                return diffCont;
            case H2OIdx:
                DUNE_THROW(Dune::InvalidStateException,
                           "Diffusivity of water in the water phase "
                           "is constraint by sum of diffusive fluxes = 0 !\n");
            }
        } else if (phaseIdx==nPhaseIdx) {
            return 0.0;
        }
        return 0;
    }

    using Base::binaryDiffusionCoefficient;
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::H2OAirXylene::binaryDiffusionCoefficient()");
    }

    using Base::fugacityCoefficient;
    /*!
     * \brief Returns the fugacity coefficient \f$\mathrm{[-]}\f$ of a component in a
     *        phase.
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase to consider
     * \param compIdx The index of the component to consider
     *
     * In this case, things are actually pretty simple. We have an ideal
     * solution. Thus, the fugacity coefficient is 1 in the gas phase
     * (fugacity equals the partial pressure of the component in the gas phase)
     * respectively in the liquid phases it is the Henry coefficients
     * divided by pressure.
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

        if (phaseIdx == wPhaseIdx) {
            if (compIdx == H2OIdx)
                return H2O::vaporPressure(T)/p;
            else if (compIdx == AirIdx)
                return BinaryCoeff::H2O_Air::henry(T)/p;
            else if (compIdx == NAPLIdx)
                return BinaryCoeff::H2O_Xylene::henry(T)/p;
        }

        // for the NAPL phase, we assume currently that nothing is
        // dissolved. this means that the affinity of the NAPL
        // component to the NAPL phase is much higher than for the
        // other components, i.e. the fugacity coefficient is much
        // smaller.
        if (phaseIdx == nPhaseIdx) {
            Scalar phiNapl = NAPL::vaporPressure(T)/p;
            if (compIdx == NAPLIdx)
                return phiNapl;
            else if (compIdx == AirIdx)
                return 1e6*phiNapl;
            else if (compIdx == H2OIdx)
                return 1e6*phiNapl;
        }

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        assert(phaseIdx == gPhaseIdx);
        return 1.0;
    }

    template <class FluidState>
    static Scalar kelvinVaporPressure(const FluidState &fluidState,
                                      const int phaseIdx,
                                      const int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::H2OAirXylene::kelvinVaporPressure()");
    }

    using Base::enthalpy;
    /*!
     * \brief Given all mole fractions in a phase, return the specific
     *        phase enthalpy \f$\mathrm{[J/kg]}\f$.
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase to consider
     *
     *  \note This system neglects the contribution of gas-molecules in the liquid phase.
     *        This contribution is probably not big. Somebody would have to find out the enthalpy of solution for this system. ...
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           int phaseIdx)
    {
        if (phaseIdx == wPhaseIdx) {
            return H2O::liquidEnthalpy(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        }
        else if (phaseIdx == nPhaseIdx) {
            return NAPL::liquidEnthalpy(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        }
        else if (phaseIdx == gPhaseIdx) {  // gas phase enthalpy depends strongly on composition
            Scalar hgc = NAPL::gasEnthalpy(fluidState.temperature(phaseIdx),
                                           fluidState.pressure(phaseIdx));
            Scalar hgw = H2O::gasEnthalpy(fluidState.temperature(phaseIdx),
                                          fluidState.pressure(phaseIdx));
            // pressure is only a dummy here (not dependent on pressure, just temperature)
            Scalar hga = Air::gasEnthalpy(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));

            Scalar result = 0;
            result += hgw * fluidState.massFraction(gPhaseIdx, H2OIdx);
            result += hga * fluidState.massFraction(gPhaseIdx, AirIdx);
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
            else if (componentIdx == AirIdx)
                DUNE_THROW(Dune::NotImplemented, "The component enthalpy for Air in water is not implemented.");
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << componentIdx);
        }
        else if (phaseIdx == nPhaseIdx)
        {
            if (componentIdx == H2OIdx)
                DUNE_THROW(Dune::NotImplemented, "The component enthalpy for water in NAPL is not implemented.");
            else if (componentIdx == NAPLIdx)
                return NAPL::liquidEnthalpy(T, p);
            else if (componentIdx == AirIdx)
                DUNE_THROW(Dune::NotImplemented, "The component enthalpy for air in NAPL is not implemented.");
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << componentIdx);
        }
        else if (phaseIdx == gPhaseIdx)
        {
            if (componentIdx == H2OIdx)
                return H2O::gasEnthalpy(T, p);
            else if (componentIdx == NAPLIdx)
                return NAPL::gasEnthalpy(T, p);
            else if (componentIdx == AirIdx)
                return Air::gasEnthalpy(T,p);
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << componentIdx);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    using Base::heatCapacity;
    /*!
     * \brief Return the heat capacity in \f$\mathrm{[J/(kg K)]}\f$.
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::H2OAirXylene::heatCapacity()");
    }

    using Base::thermalConductivity;
    /*!
     * \brief Return the thermal conductivity \f$\mathrm{[W/(m K)]}\f$.
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase
     */
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
            return NAPL::liquidThermalConductivity(temperature, pressure);
        }
        else if (phaseIdx == gPhaseIdx)
        {
            return Air::gasThermalConductivity(temperature, pressure);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

private:
    static Scalar waterPhaseDensity_(Scalar T, Scalar pw, Scalar xww, Scalar xwa, Scalar xwc)
    {
        Scalar rholH2O = H2O::liquidDensity(T, pw);
        Scalar clH2O = rholH2O/H2O::molarMass();

        // this assumes each dissolved molecule displaces exactly one
        // water molecule in the liquid
        return
            clH2O*(xww*H2O::molarMass() + xwa*Air::molarMass() + xwc*NAPL::molarMass());
    }

    static Scalar gasPhaseDensity_(Scalar T, Scalar pg, Scalar xgw, Scalar xga, Scalar xgc)
    {
        return H2O::gasDensity(T, pg*xgw) + Air::gasDensity(T, pg*xga) + NAPL::gasDensity(T, pg*xgc);
    }

    static Scalar NAPLPhaseDensity_(Scalar T, Scalar pn)
    {
        return
            NAPL::liquidDensity(T, pn);
    }

};

} // end namespace FluidSystems
} // end namespace Dumux

#endif
