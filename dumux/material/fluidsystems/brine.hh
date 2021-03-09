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
 * \brief A fluid system for brine, i.e. H2O with dissolved NaCl.
 */
#ifndef DUMUX_BRINE_FLUID_SYSTEM_HH
#define DUMUX_BRINE_FLUID_SYSTEM_HH

#include <dune/common/math.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/io/name.hh>
#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/constants.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/nacl.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

namespace Dumux {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 * \brief A compositional single phase fluid system consisting of
 *        two components, which are H2O and NaCl.
 */
template< class Scalar, class H2OType = Components::TabulatedComponent<Dumux::Components::H2O<Scalar>> >
class Brine : public Base< Scalar, Brine<Scalar, H2OType>>
{
    using ThisType = Brine<Scalar, H2OType>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;

public:
    //! export the involved components
    using H2O = H2OType;
    using NaCl = Dumux::Components::NaCl<Scalar>;

    static const int numPhases = 1;     //!< Number of phases in the fluid system
    static const int numComponents = 2; //!< Number of components in the fluid system (H2O, NaCl)

    static constexpr int phase0Idx = 0; //!< Index of the first (and only) phase
    static constexpr int liquidPhaseIdx = phase0Idx; //!< The one considered phase is liquid

    static constexpr int H2OIdx = 0;         //!< index of the water component
    static constexpr int NaClIdx = 1;        //!< index of the NaCl component
    static constexpr int comp0Idx = H2OIdx;  //!< index of the first component
    static constexpr int comp1Idx = NaClIdx; //!< index of the second component

    /*!
     * \brief Return the human readable name of the phase
     * \param phaseIdx The index of the fluid phase to consider
     */
    static const std::string phaseName(int phaseIdx = liquidPhaseIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);
        return IOName::liquidPhase();
    }

    /*!
     * \brief Returns whether the fluids are miscible
     * \note There is only one phase, so miscibility makes no sense here
     */
    static constexpr bool isMiscible()
    {
        return false;
    }

    /*!
     * \brief Return whether a phase is gaseous
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isGas(int phaseIdx = liquidPhaseIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);
        return false;
    }

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
    static bool isIdealMixture(int phaseIdx = liquidPhaseIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);
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
    static bool isCompressible(int phaseIdx = liquidPhaseIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);
        return H2O::liquidIsCompressible();
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal gas.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */

    static bool isIdealGas(int phaseIdx = liquidPhaseIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);
        return false; /*we're a liquid!*/
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    /*!
     * \brief Return the human readable name of a component
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        switch (compIdx)
        {
            case H2OIdx: return H2O::name();
            case NaClIdx: return NaCl::name();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the molar mass of a component in \f$\mathrm{[kg/mol]}\f$.
     * \param compIdx The index of the component to consider
     */
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx)
        {
            case H2OIdx: return H2O::molarMass();
            case NaClIdx: return NaCl::molarMass();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /****************************************
     * thermodynamic relations
     ****************************************/
    /*!
     * \brief Initialize the fluid system's static parameters generically
     * \note If a tabulated H2O component is used, we do our best to create
     *       tables that always work.
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
            std::cout << "Initializing tables for the H2O fluid properties ("
                      << nTemp*nPress
                      << " entries).\n";

            H2O::init(tempMin, tempMax, nTemp, pressMin, pressMax, nPress);
        }
    }

    using Base::density;
    /*!
     * \brief Return the phase density [kg/m^3].
     * \note The density is compuated as a function of the salt mass fraction
     */
    template <class FluidState>
    static Scalar density(const FluidState& fluidState, int phaseIdx = liquidPhaseIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);
        const Scalar temperature = fluidState.temperature(phaseIdx);
        const Scalar pressure = fluidState.pressure(phaseIdx);
        const Scalar xNaCl = fluidState.massFraction(phaseIdx, NaClIdx);

        using std::max;
        const Scalar TempC = temperature - 273.15;
        const Scalar pMPa = pressure/1.0E6;
        const Scalar salinity = max(0.0, xNaCl);

        const Scalar rhow = H2O::liquidDensity(temperature, pressure);
        const Scalar density = rhow + 1000*salinity*(0.668 +
                                                     0.44*salinity +
                                                     1.0E-6*(300*pMPa -
                                                             2400*pMPa*salinity +
                                                             TempC*(80.0 +
                                                                    3*TempC -
                                                                    3300*salinity -
                                                                    13*pMPa +
                                                                    47*pMPa*salinity)));
        assert(density > 0.0);
        return density;
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

    using Base::viscosity;
    /*!
     * \brief Return the viscosity of the phase.
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState& fluidState, int phaseIdx = liquidPhaseIdx)
    {
        assert(phaseIdx == liquidPhaseIdx);
        const Scalar temperature = fluidState.temperature(phaseIdx);
        const Scalar xNaCl = fluidState.massFraction(phaseIdx, NaClIdx);

        using std::pow;
        using Dune::power;
        using std::exp;
        using std::max;
        const Scalar T = max(temperature, 275.0);
        const Scalar salinity = max(0.0, xNaCl);

        const Scalar T_C = T - 273.15;
        const Scalar A = ((0.42*power((pow(salinity, 0.8)-0.17), 2)) + 0.045)*pow(T_C, 0.8);
        const Scalar mu_brine = 0.1 + (0.333*salinity) + (1.65+(91.9*salinity*salinity*salinity))*exp(-A); // [cP]
        assert(mu_brine > 0.0);
        return mu_brine/1000.0; // [Pa·s]
    }

    /*!
     * \brief Vapor pressure of a component \f$\mathrm{[Pa]}\f$.
     *
     * \param fluidState The fluid state
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar vaporPressure(const FluidState& fluidState, int compIdx)
    {
        if (compIdx == H2OIdx)
        {
            // simplified version of Eq 2.29 in Vishal Jambhekar's Promo
            const Scalar temperature = fluidState.temperature(H2OIdx);
            return H2O::vaporPressure(temperature)/fluidState.massFraction(phase0Idx, H2OIdx);
        }
        else if (compIdx == NaClIdx)
            DUNE_THROW(Dune::NotImplemented, "NaCl::vaporPressure(t)");
        else
            DUNE_THROW(Dune::NotImplemented, "Invalid component index " << compIdx);
    }

    using Base::enthalpy;
    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return its specific enthalpy \f$\mathrm{[J/kg]}\f$.
     *
     * \param fluidState The fluid state
     * \param phaseIdx The index of the phase to consider
     *
     * Equations given in:
     * - Palliser & McKibbin (1998) \cite palliser1998 <BR>
     * - Michaelides (1981) \cite michaelides1981 <BR>
     * - Daubert & Danner (1989) \cite daubert1989
     *
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState& fluidState, int phaseIdx)
    {
        /*Numerical coefficients from PALLISER*/
        static const Scalar f[] = { 2.63500E-1, 7.48368E-6, 1.44611E-6, -3.80860E-10 };

        /*Numerical coefficients from MICHAELIDES for the enthalpy of brine*/
        static const Scalar a[4][3] = { { +9633.6, -4080.0, +286.49 },
                                        { +166.58, +68.577, -4.6856 },
                                        { -0.90963, -0.36524, +0.249667E-1 },
                                        { +0.17965E-2, +0.71924E-3, -0.4900E-4 } };

        const Scalar p = fluidState.pressure(phaseIdx);
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar theta = T - 273.15;
        const Scalar salSat = f[0] + f[1]*theta + f[2]*theta*theta + f[3]*theta*theta*theta;

        /*Regularization*/
        using std::min;
        using std::max;
        const Scalar xNaCl = fluidState.massFraction(phase0Idx, NaClIdx);
        const Scalar salinity = min(max(xNaCl,0.0), salSat);

        const Scalar hw = H2O::liquidEnthalpy(T, p)/1E3; /* kJ/kg */

        /*DAUBERT and DANNER*/
        const Scalar h_NaCl = (3.6710E4*T + 0.5*(6.2770E1)*T*T - ((6.6670E-2)/3)*T*T*T
                              + ((2.8000E-5)/4)*(T*T*T*T))/(58.44E3)- 2.045698e+02; /* U [kJ/kg] */

        const Scalar m = (1E3/58.44)*(salinity/(1-salinity));

        using Dune::power;
        Scalar d_h = 0;
        for (int i = 0; i<=3; i++)
            for (int j=0; j<=2; j++)
                d_h = d_h + a[i][j] * power(theta, i) * power(m, j);

        /* heat of dissolution for halite according to Michaelides 1971 */
        const Scalar delta_h = (4.184/(1E3 + (58.44 * m)))*d_h;

        /* Enthalpy of brine without any dissolved gas */
        const Scalar h_ls1 =(1-salinity)*hw + salinity*h_NaCl + salinity*delta_h; /* kJ/kg */
        return h_ls1*1E3; /*J/kg*/
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
        const Scalar T = fluidState.temperature(liquidPhaseIdx);
        const Scalar p = fluidState.pressure(liquidPhaseIdx);

        if (phaseIdx == liquidPhaseIdx)
        {
            if (componentIdx == H2OIdx)
                return H2O::liquidEnthalpy(T, p);
            else if (componentIdx == NaClIdx)
                DUNE_THROW(Dune::NotImplemented, "The component enthalpy for NaCl is not implemented.");
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << componentIdx);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    using Base::molarDensity;
    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$ of
     *        the fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density for the simple relation is defined by the
     * mass density \f$\rho_\alpha\f$ and the molar mass of the main component \f$M_\kappa\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{M_\kappa} \;.\f]
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState& fluidState, int phaseIdx = liquidPhaseIdx)
    {
        return density(fluidState, phaseIdx)/fluidState.averageMolarMass(phaseIdx);
    }

    using Base::diffusionCoefficient;
    /*!
     * \brief Returns the diffusion coefficient \f$\mathrm{[-]}\f$
     *        of a component in a phase.
     */
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState& fluidState, int phaseIdx, int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::Brine::diffusionCoefficient()");
    }

    using Base::binaryDiffusionCoefficient;
    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for components
     *        \f$\mathrm{i}\f$ and \f$\mathrm{j}\f$ in this phase.
     * \param fluidState The fluid state
     * \param phaseIdx Index of the fluid phase
     * \param compIIdx Index of the component i
     * \param compJIdx Index of the component j
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState& fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        if (phaseIdx == liquidPhaseIdx)
        {
            if (compIIdx > compJIdx)
            {
                using std::swap;
                swap(compIIdx, compJIdx);
            }
            //! \todo TODO implement binary coefficients
            // http://webserver.dmt.upm.es/~isidoro/dat1/Mass%20diffusivity%20data.pdf
            // The link above was given as a reference in brine_air fluid system.
            if (compJIdx == NaClIdx)
                return 0.12e-9;
            else
                DUNE_THROW(Dune::NotImplemented, "Binary diffusion coefficient of components "
                                                 << compIIdx << " and " << compJIdx
                                                 << " in phase " << phaseIdx);
         }

         DUNE_THROW(Dune::InvalidStateException, "Invalid phase index: " << phaseIdx);
    }

    using Base::thermalConductivity;
    /*!
     * \brief Thermal conductivity of a fluid phase \f$\mathrm{[W/(m K)]}\f$.
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     *
     * \todo TODO: For the thermal conductivity of the phases the contribution of
     *       NaCl is neglected. This contribution is probably not big, but somebody
     *       would have to find out its influence.
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState& fluidState, int phaseIdx)
    {
        if (phaseIdx == liquidPhaseIdx)
            return H2O::liquidThermalConductivity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index: " << phaseIdx);
    }

    using Base::heatCapacity;
    /*!
     * \brief Specific isobaric heat capacity of a fluid phase.
     *        \f$\mathrm{[J/(kg*K)}\f$.
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     *
     * \todo TODO: The calculation of the isobaric heat capacity is preliminary. A better
     *       description of the influence of NaCl on the phase property has to be found.
     */
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState, int phaseIdx)
    {
        if (phaseIdx == liquidPhaseIdx)
            return H2O::liquidHeatCapacity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }
};

} // end namespace FluidSystems
} // end namespace Dumux

#endif
