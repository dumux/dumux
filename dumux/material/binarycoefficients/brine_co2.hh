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
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for CO2 and brine.
 */
#ifndef DUMUX_BINARY_COEFF_BRINE_CO2_HH
#define DUMUX_BINARY_COEFF_BRINE_CO2_HH

#include <dune/common/math.hh>

#include <dumux/common/parameters.hh>
#include <dumux/material/components/brine.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/co2.hh>
#include <dumux/material/idealgas.hh>

namespace Dumux::BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for brine and CO2.
 */
template<class Scalar, class CO2Tables, bool verbose = true>
class Brine_CO2 {
    using H2O = Dumux::Components::H2O<Scalar>;
    using CO2 = Dumux::Components::CO2<Scalar, CO2Tables>;
    using IdealGas = Dumux::IdealGas<Scalar>;
    static constexpr int lPhaseIdx = 0; // index of the liquid phase
    static constexpr int gPhaseIdx = 1; // index of the gas phase

public:
    /*!
     * \brief Binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ of water in the CO2 phase.
     *
     * According to B. Xu et al. (2002) \cite xu2003 <BR>
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        static const bool hasGasDiffCoeff = hasParam("BinaryCoefficients.GasDiffCoeff");
        if (!hasGasDiffCoeff) //in case one might set that user-specific as e.g. in dumux-lecture/mm/convectivemixing
        {
            //Diffusion coefficient of water in the CO2 phase
            constexpr Scalar PI = 3.141593;
            constexpr Scalar k = 1.3806504e-23; // Boltzmann constant
            constexpr Scalar c = 4; // slip parameter, can vary between 4 (slip condition) and 6 (stick condition)
            constexpr Scalar R_h = 1.72e-10; // hydrodynamic radius of the solute
            const Scalar mu = CO2::gasViscosity(temperature, pressure); // CO2 viscosity
            return k / (c * PI * R_h) * (temperature / mu);
        }
        else
        {
            static const Scalar D = getParam<Scalar>("BinaryCoefficients.GasDiffCoeff");
            return D;
        }
    }

    /*!
     * \brief Binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ of CO2 in the brine phase.
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        //Diffusion coefficient of CO2 in the brine phase
        static const bool hasLiquidDiffCoeff = hasParam("BinaryCoefficients.LiquidDiffCoeff");
        if (!hasLiquidDiffCoeff) //in case one might set that user-specific as e.g. in dumux-lecture/mm/convectivemixing
            return 2e-9;
        else
        {
            const Scalar D = getParam<Scalar>("BinaryCoefficients.LiquidDiffCoeff");
            return D;
        }
    }

    /*!
     * \brief Returns the _mol_ (!) fraction of CO2 in the liquid
     *        phase and the _mol_ (!) fraction of H2O in the gas phase
     *        for a given temperature, pressure, CO2 density and brine
     *        salinity.
     *
     *        Implemented according to Spycher and Pruess (2005) \cite spycher2005 <BR>
     *        applying the activity coefficient expression of Duan and Sun (2003) \cite duan2003 <BR>
     *        and the correlations for pure water given in Spycher, Pruess and Ennis-King (2003) \cite spycher2003 <BR>
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     * \param salinity the salinity \f$\mathrm{[kg \ NaCl / kg \ solution]}\f$
     * \param knownPhaseIdx indicates which phases are present
     * \param xlCO2 mole fraction of CO2 in brine \f$\mathrm{[mol/mol]}\f$
     * \param ygH2O mole fraction of water in the gas phase \f$\mathrm{[mol/mol]}\f$
     */
    static void calculateMoleFractions(const Scalar temperature,
                                       const Scalar pg,
                                       const Scalar salinity,
                                       const int knownPhaseIdx,
                                       Scalar &xlCO2,
                                       Scalar &ygH2O)
    {

        const Scalar A = computeA_(temperature, pg);

        /* salinity: conversion from mass fraction to mol fraction */
        const Scalar x_NaCl = salinityToMoleFrac_(salinity);

        // if both phases are present the mole fractions in each phase can be calculate
        // with the mutual solubility function
        if (knownPhaseIdx < 0)
        {
            const Scalar molalityNaCl = molFracToMolality_(x_NaCl); // molality of NaCl //CHANGED
            const Scalar m0_CO2 = molalityCO2inPureWater_(temperature, pg); // molality of CO2 in pure water
            const Scalar gammaStar = activityCoefficient_(temperature, pg, molalityNaCl);// activity coefficient of CO2 in brine
            const Scalar m_CO2 = m0_CO2 / gammaStar; // molality of CO2 in brine
            xlCO2 = m_CO2 / (molalityNaCl + 55.508 + m_CO2); // mole fraction of CO2 in brine
            ygH2O = A * (1 - xlCO2 - x_NaCl); // mole fraction of water in the gas phase
        }

        // if only liquid phase is present the mole fraction of CO2 in brine is given and
        // and the virtual equilibrium mole fraction of water in the non-existing gas phase can be estimated
        // with the mutual solubility function
        if (knownPhaseIdx == lPhaseIdx)
            ygH2O = A * (1 - xlCO2 - x_NaCl);

        // if only gas phase is present the mole fraction of water in the gas phase is given and
        // and the virtual equilibrium mole fraction of CO2 in the non-existing liquid phase can be estimated
        // with the mutual solubility function
        if (knownPhaseIdx == gPhaseIdx)
            xlCO2 = 1 - x_NaCl - ygH2O / A;
    }

    /*!
     * \brief Returns the fugacity coefficient of the CO2 component in a water-CO2 mixture
     * (given in Spycher, Pruess and Ennis-King (2003) \cite spycher2003 )
     *
     * \param T the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar fugacityCoefficientCO2(Scalar T, Scalar pg)
    {
        const Scalar V = 1 / (CO2::gasDensity(T, pg) / CO2::molarMass()) * 1.e6; // molar volume in cm^3/mol
        const Scalar pg_bar = pg / 1.e5; // gas phase pressure in bar
        const Scalar a_CO2 = (7.54e7 - 4.13e4 * T); // mixture parameter of  Redlich-Kwong equation
        constexpr Scalar b_CO2 = 27.8; // mixture parameter of Redlich-Kwong equation
        constexpr Scalar R = IdealGas::R * 10.; // ideal gas constant with unit bar cm^3 /(K mol)

        using std::log;
        using std::exp;
        using std::pow;
        const Scalar lnPhiCO2 = log(V / (V - b_CO2)) + b_CO2 / (V - b_CO2) - 2 * a_CO2 / (R
                * pow(T, 1.5) * b_CO2) * log((V + b_CO2) / V) + a_CO2 * b_CO2
                / (R * pow(T, 1.5) * b_CO2 * b_CO2) * (log((V + b_CO2) / V)
                - b_CO2 / (V + b_CO2)) - log(pg_bar * V / (R * T));

        return exp(lnPhiCO2); // fugacity coefficient of CO2
    }

    /*!
     * \brief Returns the fugacity coefficient of the H2O component in a water-CO2 mixture
     * (given in Spycher, Pruess and Ennis-King (2003) \cite spycher2003 )
     *
     * \param T the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar fugacityCoefficientH2O(Scalar T, Scalar pg)
    {
        const Scalar V = 1 / (CO2::gasDensity(T, pg) / CO2::molarMass()) * 1.e6; // molar volume in cm^3/mol
        const Scalar pg_bar = pg / 1.e5; // gas phase pressure in bar
        const Scalar a_CO2 = (7.54e7 - 4.13e4 * T);// mixture parameter of  Redlich-Kwong equation
        constexpr Scalar a_CO2_H2O = 7.89e7;// mixture parameter of Redlich-Kwong equation
        constexpr Scalar b_CO2 = 27.8;// mixture parameter of Redlich-Kwong equation
        constexpr Scalar b_H2O = 18.18;// mixture parameter of Redlich-Kwong equation
        constexpr Scalar R = IdealGas::R * 10.; // ideal gas constant with unit bar cm^3 /(K mol)

        using std::log;
        using std::pow;
        using std::exp;
        const Scalar lnPhiH2O = log(V / (V - b_CO2)) + b_H2O / (V - b_CO2) - 2 * a_CO2_H2O
                / (R * pow(T, 1.5) * b_CO2) * log((V + b_CO2) / V) + a_CO2
                * b_H2O / (R * pow(T, 1.5) * b_CO2 * b_CO2) * (log((V + b_CO2)
                / V) - b_CO2 / (V + b_CO2)) - log(pg_bar * V / (R * T));

        return exp(lnPhiH2O); // fugacity coefficient of H2O
    }

private:
    /*!
     * \brief Returns the molality of NaCl \f$\mathrm{[mol \ NaCl / kg \ water]}\f$  for a given mole fraction
     * \param salinity the salinity \f$\mathrm{[kg \ NaCl / kg \ solution]}\f$
     */
    static Scalar salinityToMoleFrac_(Scalar salinity)
    {
        constexpr Scalar Mw = H2O::molarMass(); // molecular weight of water [kg/mol]
        constexpr Scalar Ms = 58.8e-3;          // molecular weight of NaCl  [kg/mol]

        const Scalar X_NaCl = salinity;
        // salinity: conversion from mass fraction to mol fraction
        const Scalar x_NaCl = -Mw * X_NaCl / ((Ms - Mw) * X_NaCl - Ms);
        return x_NaCl;
    }

    /*!
     * \brief Returns the molality of NaCl \f$\mathrm{(mol \ NaCl / kg \ water)}\f$
     *        for a given mole fraction \f$\mathrm{(mol \ NaCl / mol\  solution)}\f$
     *
     * \param x_NaCl mole fraction of NaCL in brine \f$\mathrm{[mol/mol]}\f$
     */
    static Scalar molFracToMolality_(Scalar x_NaCl)
    {
        // conversion from mol fraction to molality (dissolved CO2 neglected)
        return 55.508 * x_NaCl / (1 - x_NaCl);
    }

    /*!
     * \brief Returns the equilibrium molality of CO2 \f$\mathrm{(mol \ CO2 / kg \ water)}\f$ for a
     * CO2-water mixture at a given pressure and temperature
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar molalityCO2inPureWater_(Scalar temperature, Scalar pg)
    {
        const Scalar A = computeA_(temperature, pg); // according to Spycher, Pruess and Ennis-King (2003)
        const Scalar B = computeB_(temperature, pg); // according to Spycher, Pruess and Ennis-King (2003)
        const Scalar yH2OinGas = (1 - B) / (1. / A - B); // equilibrium mol fraction of H2O in the gas phase
        const Scalar xCO2inWater = B * (1 - yH2OinGas); // equilibrium mol fraction of CO2 in the water phase
        return (xCO2inWater * 55.508) / (1 - xCO2inWater); // CO2 molality
    }

    /*!
     * \brief Returns the activity coefficient of CO2 in brine for a
     *           molal description. According to Duan and Sun (2003) \cite duan2003 <BR>
     *           given in Spycher and Pruess (2005) \cite spycher2005 <BR>
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     * \param molalityNaCl molality of NaCl \f$\mathrm{(mol \ NaCl / kg \ water)}\f$
     */
    static Scalar activityCoefficient_(Scalar temperature, Scalar pg, Scalar molalityNaCl)
    {
        const Scalar lambda = computeLambda_(temperature, pg); // lambda_{CO2-Na+}
        const Scalar xi = computeXi_(temperature, pg); // Xi_{CO2-Na+-Cl-}
        const Scalar lnGammaStar = 2 * lambda * molalityNaCl + xi * molalityNaCl
                * molalityNaCl;
        using std::exp;
        return exp(lnGammaStar); // molal activity coefficient of CO2 in brine
    }

    /*!
     * \brief Returns the paramater A for the calculation of
     * them mutual solubility in the water-CO2 system.
     * Given in Spycher, Pruess and Ennis-King (2003) \cite spycher2003 <BR>
     *
     * \param T the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar computeA_(Scalar T, Scalar pg)
    {
        const Scalar deltaP = pg / 1e5 - 1; // pressure range [bar] from p0 = 1bar to pg[bar]
        const Scalar v_av_H2O = 18.1; // average partial molar volume of H2O [cm^3/mol]
        constexpr Scalar R = IdealGas::R * 10;
        const Scalar k0_H2O = equilibriumConstantH2O_(T); // equilibrium constant for H2O at 1 bar
        const Scalar phi_H2O = fugacityCoefficientH2O(T, pg); // fugacity coefficient of H2O for the water-CO2 system
        const Scalar pg_bar = pg / 1.e5;
        using std::exp;
        return k0_H2O / (phi_H2O * pg_bar) * exp(deltaP * v_av_H2O / (R * T));
    }

    /*!
     * \brief Returns the paramater B for the calculation of
     * the mutual solubility in the water-CO2 system.
     * Given in Spycher, Pruess and Ennis-King (2003) \cite spycher2003 <BR>
     *
     * \param T the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar computeB_(Scalar T, Scalar pg)
    {
        const Scalar deltaP = pg / 1e5 - 1; // pressure range [bar] from p0 = 1bar to pg[bar]
        constexpr Scalar v_av_CO2 = 32.6; // average partial molar volume of CO2 [cm^3/mol]
        constexpr Scalar R = IdealGas::R * 10;
        const Scalar k0_CO2 = equilibriumConstantCO2_(T); // equilibrium constant for CO2 at 1 bar
        const Scalar phi_CO2 = fugacityCoefficientCO2(T, pg); // fugacity coefficient of CO2 for the water-CO2 system
        const Scalar pg_bar = pg / 1.e5;
        using std::exp;
        return phi_CO2 * pg_bar / (55.508 * k0_CO2) * exp(-(deltaP
                * v_av_CO2) / (R * T));
    }

    /*!
     * \brief Returns the parameter lambda, which is needed for the
     * calculation of the CO2 activity coefficient in the brine-CO2 system.
     * Given in Spycher and Pruess (2005) \cite spycher2005 <BR>
     * \param T the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar computeLambda_(Scalar T, Scalar pg)
    {
        constexpr Scalar c[6] = { -0.411370585, 6.07632013E-4, 97.5347708,
                -0.0237622469, 0.0170656236, 1.41335834E-5 };

        using std::log;
        const Scalar pg_bar = pg / 1.0E5; /* conversion from Pa to bar */
        return c[0] + c[1] * T + c[2] / T + c[3] * pg_bar / T + c[4] * pg_bar
                / (630.0 - T) + c[5] * T * log(pg_bar);
    }

    /*!
     * \brief Returns the parameter xi, which is needed for the
     * calculation of the CO2 activity coefficient in the brine-CO2 system.
     * Given in Spycher and Pruess (2005) \cite spycer2005 <BR>
     * \param T the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar computeXi_(Scalar T, Scalar pg)
    {
        constexpr Scalar c[4] = { 3.36389723E-4, -1.98298980E-5,
                2.12220830E-3, -5.24873303E-3 };

        Scalar pg_bar = pg / 1.0E5; /* conversion from Pa to bar */
        return c[0] + c[1] * T + c[2] * pg_bar / T + c[3] * pg_bar / (630.0 - T);
    }

    /*!
     * \brief Returns the equilibrium constant for CO2, which is needed for the
     * calculation of the mutual solubility in the water-CO2 system
     * Given in Spycher, Pruess and Ennis-King (2003) \cite spycher2003 <BR>
     * \param T the temperature \f$\mathrm{[K]}\f$
     */
    static Scalar equilibriumConstantCO2_(Scalar T)
    {
        const Scalar TinC = T - 273.15; //temperature in °C
        constexpr Scalar c[3] = { 1.189, 1.304e-2, -5.446e-5 };
        const Scalar logk0_CO2 = c[0] + c[1] * TinC + c[2] * TinC * TinC;
        using std::pow;
        return pow(10, logk0_CO2);
    }

    /*!
     * \brief Returns the equilibrium constant for H2O, which is needed for the
     * calculation of the mutual solubility in the water-CO2 system
     * Given in Spycher, Pruess and Ennis-King (2003) \cite spycher2003 <BR>
     * \param T the temperature \f$\mathrm{[K]}\f$
     */
    static Scalar equilibriumConstantH2O_(Scalar T)
    {
        const Scalar TinC = T - 273.15; //temperature in °C
        constexpr Scalar c[4] = { -2.209, 3.097e-2, -1.098e-4, 2.048e-7 };
        const Scalar logk0_H2O = c[0] + c[1] * TinC + c[2] * TinC * TinC + c[3]
                * TinC * TinC * TinC;
        using std::pow;
        return pow(10, logk0_H2O);
    }
};

/*!
 * \brief Old version of binary coefficients for CO2 and brine.
 * Calculates molfraction of CO2 in brine according to Duan and Sun 2003
 * molfraction of H2O has been assumed to be a constant value
 * For use with the actual brine_co2_system this class still needs to be adapted
 */
template<class Scalar, class CO2Tables, bool verbose = true>
class Brine_CO2_Old
{
    using H2O = Dumux::Components::H2O<Scalar>;
    using Brine = Dumux::Components::Brine<Scalar,H2O>;
    using CO2 = Dumux::Components::CO2<Scalar, CO2Tables>;
    using IdealGas = Dumux::IdealGas<Scalar>;

public:
    /*!
     * \brief Returns the _mole_ (!) fraction of CO2 in the liquid
     *        phase at a given temperature, pressure and density of
     *        CO2.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     * \param rhoCO2 density of CO2
     */
    static Scalar moleFracCO2InBrine(Scalar temperature, Scalar pg, Scalar rhoCO2)
    {
        // regularisations:
        if (pg > 2.5e8) {
            pg = 2.5e8;
        }
        if (pg < 2.e5) {
            pg = 2.e5;
        }
        if (temperature < 275.) {
            temperature = 275;
        }
        if (temperature > 600.) {
            temperature = 600;
        }

        const Scalar Mw = H2O::molarMass(); /* molecular weight of water [kg/mol] */
        const Scalar Ms = 58.8e-3; /* molecular weight of NaCl  [kg/mol] */

        const Scalar X_NaCl = Brine::salinity();
        /* salinity: conversion from mass fraction to mole fraction */
        const Scalar x_NaCl = -Mw * X_NaCl / ((Ms - Mw) * X_NaCl - Ms);

        // salinity: conversion from mole fraction to molality
        const Scalar mol_NaCl = -55.56 * x_NaCl / (x_NaCl - 1);

        const Scalar A = computeA_(temperature, pg); /* mu_{CO2}^{l(0)}/RT */
        const Scalar B = computeB_(temperature, pg); /* lambda_{CO2-Na+} */
        const Scalar C = computeC_(temperature, pg); /* Xi_{CO2-Na+-Cl-} */
        const Scalar pgCO2 = partialPressureCO2_(temperature, pg);
        const Scalar phiCO2 = fugacityCoeffCO2_(temperature, pgCO2, rhoCO2);

        using std::log;
        using Dune::power;
        const Scalar exponent = A - log(phiCO2) + 2*B*mol_NaCl + C*power(mol_NaCl,2);

        using std::exp;
        const Scalar mol_CO2w = pgCO2 / (1e5 * exp(exponent)); /* paper: equation (6) */

        const Scalar x_CO2w = mol_CO2w / (mol_CO2w + 55.56); /* conversion: molality to mole fraction */
        //Scalar X_CO2w = x_CO2w*MCO2/(x_CO2w*MCO2 + (1-x_CO2w)*Mw);   /* conversion: mole fraction to mass fraction */

        return x_CO2w;
    }

private:
    /*!
     * \brief computation of \f$\mathrm{[mu_{CO2}^{l(0)}/RT]}\f$
     * \param T the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar computeA_(Scalar T, Scalar pg)
    {
        static const Scalar c[10] = {
            28.9447706,
            -0.0354581768,
            -4770.67077,
            1.02782768E-5,
            33.8126098,
            9.04037140E-3,
            -1.14934031E-3,
            -0.307405726,
            -0.0907301486,
            9.32713393E-4,
        };

        const Scalar pg_bar = pg / 1.0E5; /* conversion from Pa to bar */
        const Scalar Tr = 630.0 - T;

        using std::log;
        return
            c[0] +
            c[1]*T +
            c[2]/T +
            c[3]*T*T +
            c[4]/Tr +
            c[5]*pg_bar +
            c[6]*pg_bar*log(T) +
            c[7]*pg_bar/T +
            c[8]*pg_bar/Tr +
            c[9]*pg_bar*pg_bar/(Tr*Tr);
    }

    /*!
     * \brief computation of B
     *
     * \param T the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar computeB_(Scalar T, Scalar pg)
    {
        const Scalar c1 = -0.411370585;
        const Scalar c2 = 6.07632013E-4;
        const Scalar c3 = 97.5347708;
        const Scalar c8 = -0.0237622469;
        const Scalar c9 = 0.0170656236;
        const Scalar c11 = 1.41335834E-5;

        const Scalar pg_bar = pg / 1.0E5; /* conversion from Pa to bar */

        using std::log;
        return
            c1 +
            c2*T +
            c3/T +
            c8*pg_bar/T +
            c9*pg_bar/(630.0-T) +
            c11*T*log(pg_bar);
    }

    /*!
     * \brief computation of C
     *
     * \param T the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar computeC_(Scalar T, Scalar pg)
    {
        const Scalar c1 = 3.36389723E-4;
        const Scalar c2 = -1.98298980E-5;
        const Scalar c8 = 2.12220830E-3;
        const Scalar c9 = -5.24873303E-3;

        const Scalar pg_bar = pg / 1.0E5; /* conversion from Pa to bar */

        return
            c1 +
            c2*T +
            c8*pg_bar/T +
            c9*pg_bar/(630.0-T);
    }

    /*!
     * \brief computation of partial pressure CO2
     *
     * We assume that the partial pressure of brine is its vapor pressure.
     * \warning: Strictly this is assumption is invalid for CO2 because the
     *           mole fraction of CO2 in brine can be considerable
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar partialPressureCO2_(Scalar temperature, Scalar pg)
    {
        return pg - Brine::vaporPressure(temperature);
    }

    /*!
     * \brief The fugacity coefficient of CO2 for a CO2-H2O mixture.
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     * \param rhoCO2 the density of CO2 for the critical volume \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar fugacityCoeffCO2_(Scalar temperature,
                                    Scalar pg,
                                    Scalar rhoCO2)
    {
        static const Scalar a[15] = {
            8.99288497E-2,
            -4.94783127E-1,
            4.77922245E-2,
            1.03808883E-2,
            -2.82516861E-2,
            9.49887563E-2,
            5.20600880E-4,
            -2.93540971E-4,
            -1.77265112E-3,
            -2.51101973E-5,
            8.93353441E-5,
            7.88998563E-5,
            -1.66727022E-2,
            1.3980,
            2.96000000E-2
        };

        // reduced temperature
        const Scalar Tr = temperature / CO2::criticalTemperature();
        // reduced pressure
        const Scalar pr = pg / CO2::criticalPressure();

        // reduced molar volume. ATTENTION: Vc is _NOT_ the critical
        // molar volume of CO2. See the reference!
        const Scalar Vc = IdealGas::R*CO2::criticalTemperature()/CO2::criticalPressure();
        const Scalar Vr =
        // molar volume of CO2 at (temperature, pg)
            CO2::molarMass() / rhoCO2
            *
                // "pseudo-critical" molar volume
                        1.0 / Vc;

        // the Z coefficient
        const Scalar Z = pr * Vr / Tr;

        const Scalar A = a[0] + a[1] / (Tr * Tr) + a[2] / (Tr * Tr * Tr);
        const Scalar B = a[3] + a[4] / (Tr * Tr) + a[5] / (Tr * Tr * Tr);
        const Scalar C = a[6] + a[7] / (Tr * Tr) + a[8] / (Tr * Tr * Tr);
        const Scalar D = a[9] + a[10] / (Tr * Tr) + a[11] / (Tr * Tr * Tr);

        using std::log;
        using std::exp;
        const Scalar lnphiCO2 =
            Z - 1 -
            log(Z) +
            A/Vr +
            B/(2*Vr*Vr) +
            C/(4*Vr*Vr*Vr*Vr) +
            D/(5*Vr*Vr*Vr*Vr*Vr)
            +
            a[12]/(2*Tr*Tr*Tr*a[14])*
            (
                a[13] + 1 -
                (  a[13] + 1 +
                   a[14]/(Vr*Vr)
                    )*exp(-a[14]/(Vr*Vr)));

        return exp(lnphiCO2);
    }

};

} // end namespace Dumux::BinaryCoeff

#endif
