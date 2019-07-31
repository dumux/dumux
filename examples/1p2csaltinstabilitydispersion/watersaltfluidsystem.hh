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
 * \brief A fluid system with one phase and two components
 *        (interstitial fluid and TRAIL, a therapeutic agent for
 *        cancer therapy).
 */
#ifndef DUMUX_WATER_SALT_FLUID_SYSTEM_HH
#define DUMUX_WATER_SALT_FLUID_SYSTEM_HH

#include <dune/common/exceptions.hh>

#include <dumux/material/fluidsystems/base.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/h2o.hh>
#include <assert.h>


namespace Dumux {

namespace FluidSystems {

/*!
 *
 * \brief A fluid system with one phase and two components
 *        (interstitial fluid and TRAIL, a therapeutic agent for
 *        cancer therapy).
 *
 * A fluid system with one phase and two components representing an
 * interstitial fluid that contains therapeutic agent (TRAIL). This is
 * used in conjunction the 1p2c model.
 */
template <class Scalar>
class WaterSalt : public Base<Scalar, WaterSalt<Scalar> >
{
    using ThisType = WaterSalt<Scalar>;
    using BaseT = Base<Scalar, ThisType>;
    using H2O = Components::TabulatedComponent<Components::H2O<Scalar>>;

public:
    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! Number of phases in the fluid system
    static constexpr int numPhases = 1;

    //! Index of the liquid phase
    static constexpr int lPhaseIdx = 1;

    /*!
     * \brief Return the hum1.2e-9an readable name of a fluid phase
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx)
    {
        static const std::string name[] =
        {
            std::string("liq")
        };
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    /*!
     * \brief Return whether a phase is liquid
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isLiquid(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return true;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are indepent on the fluid composition. This assumtion is true
     * if Henry's law and Rault's law apply. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isIdealMixture(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        // we assume Henry's and Rault's laws for the water phase and
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
    static bool isCompressible(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return true;
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! Number of components in the fluid system
    static constexpr int numComponents = 2;
    //! Index of component representing the interstitial fluid
    static constexpr int WaterIdx = 0;
    //! Index of component representing TRAIL
    static constexpr int SaltIdx = 1;

    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {
        static const std::string name[] = {
            std::string("Water"),
            std::string("Salt")
        };
        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    /*!
     * \brief Return the molar mass of a component in [kg/mol].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar molarMass(int compIdx)
    {
        static const Scalar M[] = {
            18e-3, // [kg/mol],
            58.44e-3, // [kg/mol]
        };

        assert(0 <= compIdx && compIdx < numComponents);
        return M[compIdx];
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
        if (H2O::isTabulated)
        {
            H2O::init(tempMin, tempMax, nTemp,
                      pressMin, pressMax, nPress);
        }
    }


    using BaseT::molarDensity;
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

        // assume pure water or that each gas molecule displaces exactly one
        // molecule in the liquid.
        return H2O::liquidMolarDensity(T, p);
    }

    /*!
     * \brief Return the phase density [kg/m^3].
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    using BaseT::density;
    template <class FluidState>
    static Scalar density(const FluidState &fluidState, int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        Scalar temperature = fluidState.temperature(phaseIdx);
        Scalar pressure = fluidState.pressure(phaseIdx);

        return liquidDensity1_(temperature,
               pressure,
               fluidState.massFraction(phaseIdx, SaltIdx));
    }

    /*!
     * \brief Calculate the fugacity coefficient [Pa] of an individual
     *        component in a fluid phase
     *
     * The fugacity coefficient \f$\phi_\kappa\f$ is connected to the
     * fugacity \f$f_\kappa\f$ and the component's molarity
     * \f$x_\kappa\f$ by means of the relation
     *
     * \f[ f_\kappa = \phi_\kappa * x_{\kappa} \f]
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    using BaseT::fugacityCoefficient;
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);
        return 1.0;
    }

    /*!
     * \brief Return the dynamic viscosity of a phase [Pa s].
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    using BaseT::viscosity;
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState, int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return 0.001; // [Pa*s]
    }

    /*!
     * \brief Calculate the molecular diffusion coefficient for a
     *        component in a fluid phase [mol^2 * s / (kg*m^3)]
     *
     * Molecular diffusion of a compoent \f$\kappa\f$ is caused by a
     * gradient of the chemical potential and follows the law
     *
     * \f[ J = - D \mathbf{grad} mu_\kappa \f]
     *
     * where \f$\mu_\kappa\f$ is the component's chemical potential,
     * \f$D\f$ is the diffusion coefficient and \f$J\f$ is the
     * diffusive flux. \f$mu_\kappa\f$ is connected to the component's
     * fugacity \f$f_\kappa\f$ by the relation
     *
     * \f[ \mu_\kappa = R T_\alpha \mathrm{ln} \frac{f_\kappa}{p_\alpha} \f]
     *
     * where \f$p_\alpha\f$ and \f$T_\alpha\f$ are the fluid phase'
     * pressure and temperature.
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    using BaseT::diffusionCoefficient;
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        // TODO!
        DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients");
    }

    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        return the binary diffusion coefficient for components
     *        \f$i\f$ and \f$j\f$ in this phase.
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIIdx The index of the first component to consider
     * \param compJIdx The index of the second component to consider
     */
    using BaseT::binaryDiffusionCoefficient;
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIIdx && compIIdx < numComponents);
        assert(0 <= compJIdx && compJIdx < numComponents);
        return 1.2e-9;//6.6e-6; // in [m^2/s] from Mufte modell
    }

    /*!
     * \brief Given a phase's composition, temperature, pressure and
     *        density, calculate its specific enthalpy [J/kg].
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx  for which phase to give back the heat capacity
     */
    using BaseT::enthalpy;
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState, int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        DUNE_THROW(Dune::NotImplemented, "Enthalpies");
    }

    /*!
     * \brief Thermal conductivity of a fluid phase [W/(m^2 K/m)].
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx  for which phase to give back the heat capacity
     */
    using BaseT::thermalConductivity;
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState, int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        DUNE_THROW(Dune::NotImplemented, "Thermal conductivities.");
    }

    /*!
     * \brief Specific isobaric heat capacity of a fluid phase.
     *        \f$\mathrm{[J/kg]}\f$.
     *
     * \param fluidState An abitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    using BaseT::heatCapacity;
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        DUNE_THROW(Dune::NotImplemented, "Heat capacities.");
    }

private:

    /***********************************************************************
     * water density with dissolved salt
     * rho_{b} = rho_w + contribution(salt)
    ***********************************************************************/
    static Scalar liquidDensity1_(Scalar T, Scalar pl, Scalar Xsalt)
    {
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(pl);
        Valgrind::CheckDefined(Xsalt);

        if (T < 273.15)
        {
            DUNE_THROW(NumericalProblem,
                       "Liquid density for Brine and CO2 is only "
                       "defined above 273.15K (is" << T << ")");
        }

        if (pl >= 2.5e8)
        {
            DUNE_THROW(NumericalProblem,
                       "Liquid density for Brine and CO2 is only "
                       "defined below 250MPa (is" << pl << ")");
        }

        Scalar salinity=Xsalt;
        Scalar TempC = T - 273.15;
        Scalar pMPa = pl/1.0E6;

        Scalar rhow = H2O::liquidDensity(T, pl);//1000;
            Scalar rho_brine=
                rhow +
                1000*salinity*(
                    0.668 +
                    0.44*salinity +
                    1.0E-6*(
                        300*pMPa -
                        2400*pMPa*salinity +
                        TempC*(
                            80.0 -
                            3*TempC -
                            3300*salinity -
                            13*pMPa +
                            47*pMPa*salinity)));

        return rho_brine;
    }

    /*!
        * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of pure brine.
        *
        * \param temperature temperature of component in \f$\mathrm{[K]}\f$
        * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
        *
        * Equation given in:    - Batzle & Wang (1992)
        *                         - cited by: Bachu & Adams (2002)
        *                           "Equations of State for basin geofluids"
        */
    static Scalar liquidViscosity_(Scalar temperature, Scalar pressure, Scalar Xsalt)
    {
        // regularisation
        if (temperature <= 275.)
               temperature = 275;

        Scalar T_C = temperature - 273.15;
        Scalar salinity=Xsalt;
        Scalar A = (0.42*pow((pow(salinity, 0.8)-0.17), 2) + 0.045)*pow(T_C, 0.8);
        Scalar mu_brine = 0.1 + 0.333*salinity + (1.65+91.9*salinity*salinity*salinity)*exp(-A);

        return mu_brine/1000.0; /* unit: Pa s */
    }
};


} // end namepace FluidSystems

#ifdef DUMUX_PROPERTIES_HH

/*!
 * \brief A pure single-phase fluid system.
 *
 * This is an adapter to use Dumux::WaterSaltFluidSystem<TypeTag>, as is
 * done with most other classes in Dumux and all template parameters
 * are usually defined in the property system anyhow.
 */
template<class TypeTag>
class WaterSaltFluidSystem
: public FluidSystems::WaterSalt<GetPropType<TypeTag, Properties::Scalar>>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ThisType =  WaterSaltFluidSystem<TypeTag>;

public:
    using ParameterCache = Dumux::NullParameterCache;

/*
        * \brief Given a phase's composition, temperature and pressure,
        *        return the binary diffusion coefficient for components
        *        \f$i\f$ and \f$j\f$ in this phase.
        *
        * \param fluidState An abitrary fluid state
        * \param phaseIdx The index of the fluid phase to consider
        * \param compIIdx The index of the first component to consider
        * \param compJIdx The index of the second component to consider
        */
template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             const ParameterCache &paramCache,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        Scalar diffusionCoefficient_;
        diffusionCoefficient_  = getParam<Scalar>("Problem.diffusionCoefficient");

        return diffusionCoefficient_;//1.2e-9;//6.6e-6; // in [m^2/s] from Mufte modell
    }

};

#endif // DUMUX_PROPERTIES_HH
} // end namespace Dumux

#endif // DUMUX_WATER_SALT_FLUID_SYSTEM_HH
