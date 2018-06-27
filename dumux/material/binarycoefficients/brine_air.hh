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
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for Brine and Air.
 */
#ifndef DUMUX_BINARY_COEFF_BRINE_AIR_HH
#define DUMUX_BINARY_COEFF_BRINE_AIR_HH

#include <dumux/material/components/brine.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/idealgas.hh>

namespace Dumux {
namespace BinaryCoeff {
/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for Brine and Air.
 */
template<class Scalar, class Air, bool verbose = true>
class Brine_Air
{
    using H2O = Dumux::Components::H2O<Scalar>;
    using IdealGas = Dumux::IdealGas<Scalar>;
    static const int wPhaseIdx = 0; // index of the liquid phase
    static const int nPhaseIdx = 1; // index of the gas phase

public:
    /*!
     *  \brief Binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ of water in the Air phase.
     *
     * According to B. Xu et al. (2003) \cite xu2003 <BR>
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        //Diffusion coefficient of water in the Air phase
        const Scalar Theta=1.8;
        const Scalar Daw=2.13e-5;  /* reference value */
        const Scalar pg0=1.e5;     /* reference pressure */
        const Scalar T0=273.15;    /* reference temperature */
        Scalar Dgaw;
        Dgaw=Daw*(pg0/pressure)*pow((temperature/T0),Theta);
        return Dgaw;
    }

    /*!
     * Lacking better data on water-air diffusion in liquids, we use at the
     * moment the diffusion coefficient of the air's main component nitrogen!!
     * \brief Diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular nitrogen in liquid water.
     *
     * The empirical equations for estimating the diffusion coefficient in
     * infinite solution which are presented in Reid, 1987 all show a
     * linear dependency on temperature. We thus simply scale the
     * experimentally obtained diffusion coefficient of Ferrell and
     * Himmelblau by the temperature.
     * \param temperature The temperature \f$\mathrm{[K]}\f$
     * \param pressure The pressure \f$\mathrm{[Pa]}\f$
     *
     * See:
     *
     * R. Reid et al. (1987, pp. 599) \cite reid1987 <BR>
     *
     * R. Ferrell, D. Himmelblau (1967, pp. 111-115) \cite ferrell1967
     */
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        //Diffusion coefficient of Air in the H2O phase
        const Scalar Texp = 273.15 + 25; // [K]
        const Scalar Dexp = 2.01e-9; // [m^2/s]
        return Dexp * temperature/Texp;
    }

    /*!
     * \brief Returns the _mol_ (!) fraction of Air in the liquid
     *        phase and the mol_ (!) fraction of H2O in the gas phase
     *        for a given temperature, pressure, Air density and brine
     *        XwNaCl.
     *
     *        Implemented according to Spycher and Pruess (2005) \cite spycher2005 <BR>
     *        applying the activity coefficient expression of Duan and Sun (2003) \cite duan2003 <BR>
     *        and the correlations for pure water given in Spycher, Pruess and Ennis-King (2003) \cite spycher2003 <BR>
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     * \param XwNaCl the XwNaCl \f$\mathrm{[kg NaCl / kg solution]}\f$
     * \param knownPhaseIdx indicates which phases are present
     * \param xwAir mole fraction of Air in brine \f$\mathrm{[mol/mol]}\f$
     * \param xnH2O mole fraction of water in the gas phase \f$\mathrm{[mol/mol]}\f$
     * \param xwNaCl the xwNaCl
     */
    static void calculateMoleFractions(const Scalar temperature,
                                       const Scalar pg,
                                       const Scalar XwNaCl,
                                       const int knownPhaseIdx,
                                       Scalar &xwAir,
                                       Scalar &xnH2O,
                                     Scalar &xwNaCl)
    {
        DUNE_THROW(Dune::InvalidStateException, "Function: " << "calculateMoleFractions" << " is invalid.");
//        Scalar A = computeA_(temperature, pg);
//
//        /* XwNaCl: conversion from mass fraction to mol fraction */
//        xwNaCl = massTomoleFrac_(XwNaCl);
//
//        // if both phases are present the mole fractions in each phase can be calculate
//        // with the mutual solubility function
//        if (knownPhaseIdx < 0) {
//            Scalar molalityNaCl = molFracToMolality_(xwNaCl); // molality of NaCl //CHANGED
//            Scalar m0_Air = molalityAirinPureWater_(temperature, pg); // molality of Air in pure water
//            Scalar gammaStar = activityCoefficient_(temperature, pg, molalityNaCl);// activity coefficient of Air in brine
//            Scalar m_Air = m0_Air / gammaStar; // molality of Air in brine
//            xwAir = m_Air / (molalityNaCl + 55.508 + m_Air); // mole fraction of Air in brine
//            xnH2O = A * (1 - xwAir - xwNaCl); // mole fraction of water in the gas phase
//        }
//
//        // if only liquid phase is present the mole fraction of Air in brine is given and
//        // and the virtual equilibrium mole fraction of water in the non-existing gas phase can be estimated
//        // with the mutual solubility function
//        if (knownPhaseIdx == wPhaseIdx) {
////            xnH2O = A * (1 - xwAir - xwNaCl);
//          DUNE_THROW(Dune::InvalidStateException, "phase index: " << "wPhaseIdx" << " is invalid.");
//
//        }
//
//        // if only gas phase is present the mole fraction of water in the gas phase is given and
//        // and the virtual equilibrium mole fraction of Air in the non-existing liquid phase can be estimated
//        // with the mutual solubility function
//        if (knownPhaseIdx == nPhaseIdx) {
//            //y_H2o = fluidstate.
////            xwAir = 1 - xwNaCl - xnH2O / A;
//          DUNE_THROW(Dune::InvalidStateException, "phase index: " << "nPhaseIdx" << " is invalid.");
//        }
    }

    /*!
     * \brief Returns the fugacity coefficient \f$\mathrm{[-]}\f$ of the Air component in a water-Air mixture
     * (given in Spycher, Pruess and Ennis-King (2003) \cite spycher2003 )
     *
     * \param T the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar fugacityCoefficientAir(Scalar T, Scalar pg)
    {
        Scalar V = 1 / (Air::gasDensity(T, pg) / Air::molarMass()) * 1.e6; // molar volume in cm^3/mol
        Scalar pg_bar = pg / 1.e5; // gas phase pressure in bar
        Scalar a_Air = (7.54e7 - 4.13e4 * T); // mixture parameter of  Redlich-Kwong equation
        static const Scalar b_Air = 27.8; // mixture parameter of Redlich-Kwong equation
        static const Scalar R = IdealGas::R * 10.; // ideal gas constant with unit bar cm^3 /(K mol)
        Scalar lnPhiAir, phiAir;

        using std::log;
        using std::pow;
        using std::exp;
        lnPhiAir = log(V / (V - b_Air)) + b_Air / (V - b_Air) - 2 * a_Air / (R
                * pow(T, 1.5) * b_Air) * log((V + b_Air) / V) + a_Air * b_Air
                / (R * pow(T, 1.5) * b_Air * b_Air) * (log((V + b_Air) / V)
                - b_Air / (V + b_Air)) - log(pg_bar * V / (R * T));

        phiAir = exp(lnPhiAir); // fugacity coefficient of Air
        return phiAir;
    }

    /*!
     * \brief Returns the fugacity coefficient \f$\mathrm{[-]}\f$ of the H2O component in a water-Air mixture
     * (given in Spycher, Pruess and Ennis-King (2003) \cite spycher2003 )
     *
     * \param T the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar fugacityCoefficientH2O(Scalar T, Scalar pg)
    {
        Scalar V = 1 / (Air::gasDensity(T, pg) / Air::molarMass()) * 1.e6; // molar volume in cm^3/mol
        Scalar pg_bar = pg / 1.e5; // gas phase pressure in bar
        Scalar a_Air = (7.54e7 - 4.13e4 * T);// mixture parameter of  Redlich-Kwong equation
        static const Scalar a_Air_H2O = 7.89e7;// mixture parameter of Redlich-Kwong equation
        static const Scalar b_Air = 27.8;// mixture parameter of Redlich-Kwong equation
        static const Scalar b_H2O = 18.18;// mixture parameter of Redlich-Kwong equation
        static const Scalar R = IdealGas::R * 10.; // ideal gas constant with unit bar cm^3 /(K mol)
        Scalar lnPhiH2O, phiH2O;

        using std::log;
        using std::pow;
        using std::exp;
        lnPhiH2O = log(V / (V - b_Air)) + b_H2O / (V - b_Air) - 2 * a_Air_H2O
                / (R * pow(T, 1.5) * b_Air) * log((V + b_Air) / V) + a_Air
                * b_H2O / (R * pow(T, 1.5) * b_Air * b_Air) * (log((V + b_Air)
                / V) - b_Air / (V + b_Air)) - log(pg_bar * V / (R * T));
        phiH2O = exp(lnPhiH2O); // fugacity coefficient of H2O
        return phiH2O;
    }

    /*!
     * \brief Returns the molality of NaCl \f$\mathrm{(mol NaCl / kg water)}\f$ for a given mole fraction \f$\mathrm{(mol NaCl / mol solution)}\f$
     *
     * \param XwNaCl mole fraction of NaCL in brine \f$\mathrm{[mol/mol]}\f$
     */
    static Scalar molalityNaCl(Scalar XwNaCl)
    {

        // conversion from mol fraction to molality
        const Scalar mol_NaCl = XwNaCl / 58.4428e-3;

        return mol_NaCl;
    }

private:
    /*!
     * \brief Returns the molality of NaCl \f$\mathrm{(mol NaCl / kg water)}\f$ for a given mole fraction
     *
     * \param XwNaCl the XwNaCl \f$\mathrm{[kg NaCl / kg solution]}\f$
     */
    static Scalar massTomoleFrac_(Scalar XwNaCl)
    {
        DUNE_THROW(Dune::InvalidStateException, "Function: " << "massTomoleFrac_" << " is invalid.");

//        const Scalar Mw = H2O::molarMass(); /* molecular weight of water [kg/mol] */
//        const Scalar Ms = 58.8e-3; /* molecular weight of NaCl  [kg/mol] */
//
//        const Scalar X_NaCl = XwNaCl;
//        /* XwNaCl: conversion from mass fraction to mol fraction */
//        const Scalar xwNaCl = -Mw * X_NaCl / ((Ms - Mw) * X_NaCl - Ms);
//        return xwNaCl;
    }

    /*!
     * \brief Returns the equilibrium molality of Air \f$\mathrm{(mol Air / kg water)}\f$ for a
     * Air-water mixture at a given pressure and temperature
     *
     * \param T the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar molalityAirinPureWater_(Scalar temperature, Scalar pg)
    {
        Scalar A = computeA_(temperature, pg); // according to Spycher, Pruess and Ennis-King (2003)
        Scalar B = computeB_(temperature, pg); // according to Spycher, Pruess and Ennis-King (2003)
        Scalar yH2OinGas = (1 - B) / (1. / A - B); // equilibrium mol fraction of H2O in the gas phase
        Scalar xAirinWater = B * (1 - yH2OinGas); // equilibrium mol fraction of Air in the water phase
        Scalar molalityAir = (xAirinWater * 55.508) / (1 - xAirinWater); // Air molality
        return molalityAir;
    }

    /*!
     * \brief Returns the activity coefficient of Air in brine for a
     *           molal description. According to Duan and Sun (2003) \cite duan2003 <BR>
     *           given in Spycher and Pruess (2005) \cite spycher2005 <BR>
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     * \param molalityNaCl molality of NaCl \f$\mathrm{(mol NaCl / kg water)}\f$
     */
    static Scalar activityCoefficient_(Scalar temperature, Scalar pg, Scalar molalityNaCl)
    {
        Scalar lambda = computeLambda_(temperature, pg); // lambda_{Air-Na+}
        Scalar xi = computeXi_(temperature, pg); // Xi_{Air-Na+-Cl-}
        Scalar lnGammaStar = 2 * lambda * molalityNaCl + xi * molalityNaCl
                * molalityNaCl;
        using std::exp;
        Scalar gammaStar = exp(lnGammaStar);
        return gammaStar; // molal activity coefficient of Air in brine
    }

    /*!
     * \brief Returns the paramater A for the calculation of
     * them mutual solubility in the water-Air system.
     * Given in Spycher, Pruess and Ennis-King (2003) \cite spycher2003 <BR>
     *
     * \param T the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar computeA_(Scalar T, Scalar pg)
    {
        Scalar deltaP = pg / 1e5 - 1; // pressure range [bar] from p0 = 1bar to pg[bar]
        const Scalar v_av_H2O = 18.1; // average partial molar volume of H2O [cm^3/mol]
        const Scalar R = IdealGas::R * 10;
        Scalar k0_H2O = equilibriumConstantH2O_(T); // equilibrium constant for H2O at 1 bar
        Scalar phi_H2O = fugacityCoefficientH2O(T, pg); // fugacity coefficient of H2O for the water-Air system
        Scalar pg_bar = pg / 1.e5;
        using std::exp;
        Scalar A = k0_H2O / (phi_H2O * pg_bar) * exp(deltaP * v_av_H2O / (R * T));
        return A;
    }

    /*!
     * \brief Returns the paramater B for the calculation of
     * the mutual solubility in the water-Air system.
     * Given in Spycher, Pruess and Ennis-King (2003) \cite spycher2003 <BR>
     *
     * \param T the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar computeB_(Scalar T, Scalar pg)
    {
        Scalar deltaP = pg / 1e5 - 1; // pressure range [bar] from p0 = 1bar to pg[bar]
        const Scalar v_av_Air = 32.6; // average partial molar volume of Air [cm^3/mol]
        const Scalar R = IdealGas::R * 10;
        Scalar k0_Air = equilibriumConstantAir_(T); // equilibrium constant for Air at 1 bar
        Scalar phi_Air = fugacityCoefficientAir(T, pg); // fugacity coefficient of Air for the water-Air system
        Scalar pg_bar = pg / 1.e5;
        using std::exp;
        Scalar B = phi_Air * pg_bar / (55.508 * k0_Air) * exp(-(deltaP
                * v_av_Air) / (R * T));
        return B;
    }

    /*!
     * \brief Returns the parameter lambda, which is needed for the
     * calculation of the Air activity coefficient in the brine-Air system.
     * Given in Spycher and Pruess (2005) \cite spycher2005 <BR>
     * \param T the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar computeLambda_(Scalar T, Scalar pg)
    {
        Scalar lambda;
        static const Scalar c[6] = { -0.411370585, 6.07632013E-4, 97.5347708,
                -0.0237622469, 0.0170656236, 1.41335834E-5 };

        using std::log;
        Scalar pg_bar = pg / 1.0E5; /* conversion from Pa to bar */
        lambda = c[0] + c[1] * T + c[2] / T + c[3] * pg_bar / T + c[4] * pg_bar
                / (630.0 - T) + c[5] * T * log(pg_bar);

        return lambda;
    }

    /*!
     * \brief Returns the parameter xi, which is needed for the
     * calculation of the Air activity coefficient in the brine-Air system.
     * Given in Spycher and Pruess (2005) \cite spycher2005 <BR>
     * \param T the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar computeXi_(Scalar T, Scalar pg)
    {
        Scalar xi;
        static const Scalar c[4] = { 3.36389723E-4, -1.98298980E-5,
                2.12220830E-3, -5.24873303E-3 };

        Scalar pg_bar = pg / 1.0E5; /* conversion from Pa to bar */
        xi = c[0] + c[1] * T + c[2] * pg_bar / T + c[3] * pg_bar / (630.0 - T);

        return xi;
    }

    /*!
     * \brief Returns the equilibrium constant for Air, which is needed for the
     * calculation of the mutual solubility in the water-Air system
     * Given in Spycher, Pruess and Ennis-King (2003) \cite spycher2003 <BR>
     * \param T the temperature \f$\mathrm{[K]}\f$
     */
    static Scalar equilibriumConstantAir_(Scalar T)
    {
        Scalar TinC = T - 273.15; //temperature in °C
        static const Scalar c[3] = { 1.189, 1.304e-2, -5.446e-5 };
        Scalar logk0_Air = c[0] + c[1] * TinC + c[2] * TinC * TinC;
        Scalar k0_Air = pow(10, logk0_Air);
        return k0_Air;
    }

    /*!
     * \brief Returns the equilibrium constant for H2O, which is needed for the
     * calculation of the mutual solubility in the water-Air system
     * Given in Spycher, Pruess and Ennis-King (2003) \cite spycher2003 <BR>
     * \param T the temperature \f$\mathrm{[K]}\f$
     */
    static Scalar equilibriumConstantH2O_(Scalar T)
    {
        Scalar TinC = T - 273.15; //temperature in °C
        static const Scalar c[4] = { -2.209, 3.097e-2, -1.098e-4, 2.048e-7 };
        Scalar logk0_H2O = c[0] + c[1] * TinC + c[2] * TinC * TinC + c[3]
                * TinC * TinC * TinC;
        Scalar k0_H2O = pow(10, logk0_H2O);
        return k0_H2O;
    }
};

/*!
 * \brief Old version of binary coefficients for Air and brine.
 * Calculates molfraction of Air in brine according to Duan and Sun 2003
 * molfraction of H2O has been assumed to be a constant value
 * For use with the actual brine_co2_system this class still needs to be adapted
 */
template<class Scalar, class Air, bool verbose = true>
class Brine_Air_Old
{
    using H2O = Dumux::Components::H2O<Scalar>;
    using Brine = Dumux::Components::Brine<Scalar,H2O>;
   // using Air = Dumux::Components::Air<Scalar>;
    using IdealGas = Dumux::IdealGas<Scalar>;

public:
    /*!
     * \brief Returns the _mole_ (!) fraction of Air in the liquid
     *        phase at a given temperature, pressure and density of
     *        Air.
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pg the gas phase pressure \f$\mathrm{[Pa]}\f$
     * \param rhoAir density of Air
     */
    static Scalar moleFracAirInBrine(Scalar temperature, Scalar pg, Scalar rhoAir)
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

        const Scalar X_NaCl = Brine::XwNaCl;
        /* XwNaCl: conversion from mass fraction to mole fraction */
        const Scalar xwNaCl = -Mw * X_NaCl / ((Ms - Mw) * X_NaCl - Ms);

        // XwNaCl: conversion from mole fraction to molality
        const Scalar mol_NaCl = -55.56 * xwNaCl / (xwNaCl - 1);

        const Scalar A = computeA_(temperature, pg); /* mu_{Air}^{l(0)}/RT */
        const Scalar B = computeB_(temperature, pg); /* lambda_{Air-Na+} */
        const Scalar C = computeC_(temperature, pg); /* Xi_{Air-Na+-Cl-} */
        const Scalar pgAir = partialPressureAir_(temperature, pg);
        const Scalar phiAir = fugacityCoeffAir_(temperature, pgAir, rhoAir);

        using std::log;
        using std::pow;
        using std::exp;
        const Scalar exponent = A - log(phiAir) + 2*B*mol_NaCl + C*pow(mol_NaCl,2);

        const Scalar mol_Airw = pgAir / (1e5 * exp(exponent)); /* paper: equation (6) */

        const Scalar x_Airw = mol_Airw / (mol_Airw + 55.56); /* conversion: molality to mole fraction */
        //Scalar X_Airw = x_Airw*MAir/(x_Airw*MAir + (1-x_Airw)*Mw);   /* conversion: mole fraction to mass fraction */

        return x_Airw;
    }

private:
    /*!
     * \brief computation of mu_{Air}^{l(0)}/RT
     * \param T the temperature [K]
     * \param pg the gas phase pressure [Pa]
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
     * \param T the temperature [K]
     * \param pg the gas phase pressure [Pa]
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
     * \param T the temperature [K]
     * \param pg the gas phase pressure [Pa]
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
     * \brief computation of partial pressure Air
     *
     * We assume that the partial pressure of brine is its vapor pressure.
     * \warning: Strictly this is assumption is invalid for air because the
     *           mole fraction of air in brine can be considerable
     *
     * \param temperature the temperature [K]
     * \param pg the gas phase pressure [Pa]
     */
    static Scalar partialPressureAir_(Scalar temperature, Scalar pg)
    {
        return pg - Brine::vaporPressure(temperature);
    }

    /*!
     * \brief The fugacity coefficient of Air for a Air-H2O mixture.
     *
     * \param temperature the temperature [K]
     * \param pg the gas phase pressure [Pa]
     * \param rhoAir the density of Air for the critical volume [kg/m^3]
     */
    static Scalar fugacityCoeffAir_(Scalar temperature,
                                    Scalar pg,
                                    Scalar rhoAir)
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
        const Scalar Tr = temperature / Air::criticalTemperature();
        // reduced pressure
        const Scalar pr = pg / Air::criticalPressure();

        // reduced molar volume. ATTENTION: Vc is _NOT_ the critical
        // molar volume of Air. See the reference!
        const Scalar Vc = IdealGas::R*Air::criticalTemperature()/Air::criticalPressure();
        const Scalar Vr =
        // molar volume of Air at (temperature, pg)
            Air::molarMass() / rhoAir
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
        const Scalar lnphiAir =
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

        return exp(lnphiAir);
    }

};
}
} // end namespace

#endif
