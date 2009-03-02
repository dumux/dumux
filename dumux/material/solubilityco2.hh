// $Id$

#ifndef DUNE_SOLUBILITYCO2_HH
#define DUNE_SOLUBILITYCO2_HH

#include"dumux/material/constrel/constrelco2.hh"

namespace Dune
{

/** \todo Please doc me! */

// from MUFTE
class SolubilityCO2
{
public:

    double SolCO2inWater(double Temp, double pg, double X_NaCl)
    {
        // regularisations:

        if(pg > 2.5e8)
        {
            pg = 2.5e8;
        }

        if(pg < 2.e5)
        {
            pg = 2.e5;
        }

        if(Temp < 275.)
        {
            Temp = 275;
        }

        if(Temp > 600.)
        {
            Temp = 600;
        }
        /* X_NaCl: salinity: mass fraction [-] */
        double x_NaCl;      /* salinity: mole fraction [-] */
        double mol_NaCl;    /* salinity: molality [mol/kg water] */
        double Mw, Ms, MCO2;
        double xCO2g;       /* mole fraction of CO2 in the gasphase */
        double X_CO2w;      /* mass fraction of CO2 in waterphase [-] */
        double x_CO2w;      /* mole fraction of CO2 in waterphase [-] */
        double mol_CO2w;    /* molality of CO2 in the waterphase [mol/kg water] */
        double pg_bar;      /* pressure in bar */
        double phiCO2;
        double A, B, C, exponent;

        Mw   = 0.018;         /* molecular weight of water [kg/mol] */
        Ms   = 0.0588;        /* molecular weight of NaCl  [kg/mol] */
        MCO2 = 0.044;        /* molecular weight of CO2  [kg/mol] */

        x_NaCl   = -Mw*X_NaCl/((Ms-Mw)*X_NaCl - Ms); /* salinity: conversion from mass fraction to mole fraction */
        mol_NaCl = -55.56*x_NaCl/(x_NaCl-1);         /* salinity: conversion from mole fraction to molality */

        pg_bar = pg/1.0E5;
        xCO2g  = xCO2inGasphase(Temp, pg);

        A = ComputeA(Temp,pg);    /* mu_{CO2}^{l(0)}/RT */
        B = ComputeB(Temp,pg);    /* lambda_{CO2-Na+} */
        C = ComputeC(Temp,pg);    /* Xi_{CO2-Na+-Cl-} */
        phiCO2 = FugacityCoeffCO2(Temp,pg);

        exponent = A - log(phiCO2) + 2*B*mol_NaCl + C*pow(mol_NaCl,2);

        mol_CO2w = xCO2g * pg_bar / exp(exponent);  /* paper: equation (6) */

        x_CO2w = mol_CO2w/(mol_CO2w + 55.56);              /* conversion: molality to mole fraction */
        X_CO2w = x_CO2w*MCO2/(x_CO2w*MCO2 + (1-x_CO2w)*Mw);   /* conversion: mole fraction to mass fraction */

        return(X_CO2w);
    }

    /********************************************************************/
    /*                                                                  */
    /*      Solubility of CO2 in brine:                                 */
    /*      Calculation of Henry coefficient after                      */
    /*      Bando et al. (2003) (J. Chem. Eng. Data 2003, 48, 576-579)  */
    /*                                                                  */
    /*      validity: pressure:    10 - 20 MPa                          */
    /*                temperature: 30 - 60°C                            */
    /*                salinity:    0.01-0.03  (mass fraction, Xsw)      */
    /*                                                                  */
    /********************************************************************/

    double Henry_Bando(double Temp, double p, double Xsw)
    {
        double p_MPa;
        double Henry_MPa, Henry_Pa;

        /* caution: pressures are needed here in [MPa], temperatures in [°C] */
        /*          salinities in mass fractions (Xsw) [-]                     */

        p_MPa  = p/1.0E6;

        Henry_MPa = 36.1*p_MPa + 3.87*Temp - 1097.1 + (196*p_MPa + 26.9*Temp - 8810.0)*Xsw;

        /* we need Henry in [Pa] */
        Henry_Pa = Henry_MPa * 1.0E6;

        return(Henry_Pa);
    }

    /****************************************************************/
    /*                                                              */
    /*      solubility of CO2 in brine                              */
    /*      calculation of Henry coefficient                        */
    /*                                                              */
    /*      after: Battistelli, Calore and Pruess                   */
    /*      (Geothermics, Vol. 26, No.4, pp. 437-464, 1997)         */
    /*                                                              */
    /****************************************************************/

    double Henry_coeff(double Temp, double molality)
    {
        int i;
        double exponent;
        double Henry;
        double kb, Kh;
        double B[6], C[5];

        /* Regularisierung */
        if (Temp > 623.15) Temp = 623.15;
        if (Temp < 273.15) Temp = 273.15;

        Temp = Temp - 273.15;   /* Temp in [°C] */


        /* calculation of kb (salting out coefficient) */

        C[0] =  1.19784E-1;
        C[1] = -7.17823E-4;
        C[2] =  4.93854E-6;
        C[3] = -1.03826E-8;
        C[4] =  1.08233E-11;

        kb = 0.0;

        for (i = 0; i < 5; i++)
        {
            kb += C[i] * pow(Temp, i);
        }

        /* calculation of Henry's constant of pure water (Kh) */

        B[0] =  7.83666E7;
        B[1] =  1.96025E6;
        B[2] =  8.20574E4;
        B[3] = -7.40674E2;
        B[4] =  2.1838;
        B[5] = -2.20999E-3;

        Kh = 0.0;

        for (i = 0; i < 6; i++)
        {
            Kh += B[i] * pow(Temp, i);
        }

        exponent = molality * kb;

        Henry = Kh * pow(10, exponent);

        return(Henry);
    }

    /****************************************************************/
    /*                                                              */
    /* Heat of solution CO2 in water                                */
    /* Duan & Sun (2003)                                            */
    /*                                                              */
    /****************************************************************/

    double HeatOfSolutionCO2inWater (double T)
    {
        double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10;
        double pw, v;
        double dh_mol, dh;

        c1 = 28.9447706;
        c2 = -0.0354581768;
        c3 = -4770.67077;
        c4 = 1.02782768E-5;
        c5 = 33.8126098;
        c6 = 9.04037140E-3;
        c7 = -1.14934031E-3;
        c8 = -0.307405726;
        c9 = -0.0907301486;
        c10= 9.32713393E-4;

        pw = pure_water_pressure (T)/1.0E5;  /* water pressure in bar */
        v  = 630.0 - T;

        /* dh_mol : heat of solution in [kJ/mol] */
        dh_mol = -0.008314467*T*T*(c2 - c3/(T*T) + 2*c4*T + c5/(v*v) + c7*pw/T - c8*pw/(T*T) + c9*pw/(v*v) + 2*c10*pw*pw/(v*v*v));

        /* dh : heat of solution in [kJ/kg] */
        dh = dh_mol/0.044;

        return dh;
    }

    /***********************************************************************/
    /*                                                                      */
    /* This function computes the mole fraction of CO2 in the gaseous       */
    /* phase using Dalton's Law (Sum of partial pressures = total pressure) */
    /*                                                                      */
    /* The partial pressure of the waterphase is assumed to be the         */
    /* saturation vapor pressure (pwsat)                                    */
    /*                                                                      */
    /***********************************************************************/

    double xCO2inGasphase(double Temp, double pg)
    {
        double pH2O;    /* partial pressure of water in the gasphase */
        double xCO2g;   /* mole fraction of CO2 in the gasphase */

        pH2O = pure_water_pressure (Temp);

        xCO2g = (pg - pH2O)/pg;

        return(xCO2g);
    }
    /**********************************************************/
    /*                                                        */
    /* ComputeA: computation of mu_{CO2}^{l(0)}/RT            */
    /*                                                        */
    /**********************************************************/

    double ComputeA(double Temp, double pg)
    {
        double T, pg_bar, p;
        double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10;
        double A;

        T      = Temp;
        pg_bar = pg/1.0E5;   /* conversion from Pa to bar */
        p      = pg_bar;

        c1 = 28.9447706;
        c2 = -0.0354581768;
        c3 = -4770.67077;
        c4 = 1.02782768E-5;
        c5 = 33.8126098;
        c6 = 9.04037140E-3;
        c7 = -1.14934031E-3;
        c8 = -0.307405726;
        c9 = -0.0907301486;
        c10= 9.32713393E-4;

        A = c1+c2*T+c3/T+c4*T*T+c5/(630.0-T)+c6*p+c7*p*log(T)+c8*p/T+c9*p/(630.0-T)+c10*p*p/(pow((630.0-T),2));

        return(A);
    }

    double iterateVr(double T, double p)
    {
        double A, B, C, D;
        double Tcrit, pcrit, Tr, pr;
        double Vr;
        double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15;
        double Vc, Vd, Ve, c, d, e, i;
        double V1,V2, e1, ps;
        int n, n1, n2;

        a1  =  8.99288497E-2;
        a2  = -4.94783127E-1;
        a3  =  4.77922245E-2;
        a4  =  1.03808883E-2;
        a5  = -2.82516861E-2;
        a6  =  9.49887563E-2;
        a7  =  5.20600880E-4;
        a8  = -2.93540971E-4;
        a9  = -1.77265112E-3;
        a10 = -2.51101973E-5;
        a11 =  8.93353441E-5;
        a12 =  7.88998563E-5;
        a13 = -1.66727022E-2;
        a14 =  1.3980;
        a15 =  2.96000000E-2;

        pcrit = 7.3825E6;
        Tcrit = 304.2;

        Tr = T/Tcrit;    /* reduced temperature */
        pr = p/pcrit;    /* reduced pressure */

        A = a1 + a2/(Tr*Tr) + a3/(Tr*Tr*Tr);
        B = a4 + a5/(Tr*Tr) + a6/(Tr*Tr*Tr);
        C = a7 + a8/(Tr*Tr) + a9/(Tr*Tr*Tr);
        D = a10 + a11/(Tr*Tr) + a12/(Tr*Tr*Tr);

        V2=V1=Vr=1e5;n1=n2=0;ps=0;
        {
            Vc = 0.1; Vd = Vc+0.0001;i=-1;
            for (n=1; n<=20; n++)
            {
                c = Tr/(pr*Vc) * (1 + A/Vc + B/(Vc*Vc) + C/(Vc*Vc*Vc*Vc) + D/(Vc*Vc*Vc*Vc*Vc) + a13/(Tr*Tr*Tr*Vc*Vc)*(a14+a15/(Vc*Vc))*exp(-a15/(Vc*Vc))) - 1.0;
                d = Tr/(pr*Vd) * (1 + A/Vd + B/(Vd*Vd) + C/(Vd*Vd*Vd*Vd) + D/(Vd*Vd*Vd*Vd*Vd) + a13/(Tr*Tr*Tr*Vd*Vd)*(a14+a15/(Vd*Vd))*exp(-a15/(Vd*Vd))) - 1.0;
                Ve = Vc - (c*(Vd-Vc))/(d-c);
                e = Tr/(pr*Ve) * (1 + A/Ve + B/(Ve*Ve) + C/(Ve*Ve*Ve*Ve) + D/(Ve*Ve*Ve*Ve*Ve) + a13/(Tr*Tr*Tr*Ve*Ve)*(a14+a15/(Ve*Ve))*exp(-a15/(Ve*Ve))) - 1.0;
                if(fabs(e)<1.0E-4) break;
                Vc=Ve;
                Vd = Vc+0.0001*Vc;
            }
            V1=Ve;Vr=V1;n1=n;e1=e;
        }
        ConstrelCO2 constrelCO2;
        ps = constrelCO2.satpressure(T);
        Vc = 0.1; Vd = 85.0; i=-1;
        for (n=1; n<=30; n++)
        {
            c = Tr/(pr*Vc) * (1 + A/Vc + B/(Vc*Vc) + C/(Vc*Vc*Vc*Vc) + D/(Vc*Vc*Vc*Vc*Vc) + a13/(Tr*Tr*Tr*Vc*Vc)*(a14+a15/(Vc*Vc))*exp(-a15/(Vc*Vc))) - 1.0;
            d = Tr/(pr*Vd) * (1 + A/Vd + B/(Vd*Vd) + C/(Vd*Vd*Vd*Vd) + D/(Vd*Vd*Vd*Vd*Vd) + a13/(Tr*Tr*Tr*Vd*Vd)*(a14+a15/(Vd*Vd))*exp(-a15/(Vd*Vd))) - 1.0;
            if (i<n) Ve = (Vc+Vd)/2;
            else Ve = Vd-((Vd-Vc)*d)/(d-c);
            e = Tr/(pr*Ve) * (1 + A/Ve + B/(Ve*Ve) + C/(Ve*Ve*Ve*Ve) + D/(Ve*Ve*Ve*Ve*Ve) + a13/(Tr*Tr*Tr*Ve*Ve)*(a14+a15/(Ve*Ve))*exp(-a15/(Ve*Ve))) - 1.0;
            if(fabs(e)<1.0E-4) break;
            if ((c*e)<0) Vd=Ve;
            if ((d*e)<0) Vc=Ve;
            i=-1.2*i;
        }

        V2=Ve;n2=n;
        if (T<=304.0 && ps>p) Vr=V2;
        else if (T>304.0) Vr=V2;
        /*    printf("%5e %5e %5e %5e %5e\n", p,test1,test2,Vtest1,Vtest2);*/
        if (e > 1.0E-4)
        {
            printf("Vorsicht: keine Loesung fuer Vr gefunden!!! \t p = %.2f bar T = %.2f K \n", p/1.0E5, T);
        }
        return(Vr);
    }

    /**********************************************************/
    /*                                                        */
    /* ComputeB: computation of lambda_{CO2-Na}               */
    /*                                                        */
    /**********************************************************/


    double ComputeB(double Temp, double pg)
    {
        double c1,c2,c3,c8,c9,c11;
        double T,pg_bar,p;
        double B;

        c1 = -0.411370585;
        c2 = 6.07632013E-4;
        c3 = 97.5347708;
        c8 = -0.0237622469;
        c9 = 0.0170656236;
        c11= 1.41335834E-5;

        T      = Temp;
        pg_bar = pg/1.0E5;   /* conversion from Pa to bar */
        p      = pg_bar;

        B = c1+c2*T+c3/T+c8*p/T+c9*p/(630.0-T)+c11*T*log(p);

        return(B);
    }
    double ComputeC(double Temp, double pg)
    {
        double c1,c2,c8,c9;
        double T,pg_bar,p;
        double C;

        c1 = 3.36389723E-4;
        c2 = -1.98298980E-5;
        c8 = 2.12220830E-3;
        c9 = -5.24873303E-3;

        T      = Temp;
        pg_bar = pg/1.0E5;   /* conversion from Pa to bar */
        p      = pg_bar;

        C = c1+c2*T+c8*p/T+c9*p/(630.0-T);

        return(C);
    }
    /************************************************************************/
    /*                                                                      */
    /* Computation of real gas factor Z as suggested in Duan & Sun (1992)   */
    /* in Appendix A                                                        */
    /*                                                                      */
    /************************************************************************/


    double ComputeZ(double Temp, double pg)
    {
        double pr, pcrit, Tr, Tcrit, Vr, R, Z;

        pcrit = 7.3825E6;
        Tcrit = 304.2;
        R   = 8.314467;   /* universal gas constant [Pa m^3/(K mol)] */

        pr = pg/pcrit;    /* computation of reduced pressure */
        Tr = Temp/Tcrit;  /* computation of reduced temperature */

        Vr = iterateVr(Temp, pg);

        Z = pr*Vr/Tr;

        /*   printf("p = %.2f bar T = %.2f C \t Vr = %.4f \t Z = %.4f \n", pg/1.0E5, Temp-273.15, Vr, Z);*/


        return Z;
    }
    double FugacityCoeffCO2(double Temp, double pg)
    {
        double A, B, C, D;
        double Tcrit, Tr, pcrit, pr, Vr;
        double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15;
        double Z;
        double lnphiCO2, phiCO2;

        Z = ComputeZ(Temp, pg);

        a1  =  8.99288497E-2;
        a2  = -4.94783127E-1;
        a3  =  4.77922245E-2;
        a4  =  1.03808883E-2;
        a5  = -2.82516861E-2;
        a6  =  9.49887563E-2;
        a7  =  5.20600880E-4;
        a8  = -2.93540971E-4;
        a9  = -1.77265112E-3;
        a10 = -2.51101973E-5;
        a11 =  8.93353441E-5;
        a12 =  7.88998563E-5;
        a13 = -1.66727022E-2;
        a14 =  1.3980;
        a15 =  2.96000000E-2;

        pcrit = 7.3825E6;     /* critical pressure [Pa] */
        Tcrit = 304.2;     /* critical temperature [K] */

        Tr = Temp/Tcrit;    /* reduced temperature */
        pr = pg/pcrit;    /* reduced pressure */

        Vr = Z*Tr/pr;    /* compute Vr backwards */

        A = a1 + a2/(Tr*Tr) + a3/(Tr*Tr*Tr);
        B = a4 + a5/(Tr*Tr) + a6/(Tr*Tr*Tr);
        C = a7 + a8/(Tr*Tr) + a9/(Tr*Tr*Tr);
        D = a10 + a11/(Tr*Tr) + a12/(Tr*Tr*Tr);

        lnphiCO2 =  Z - 1 - log(Z) + A/Vr + B/(2*Vr*Vr) + C/(4*Vr*Vr*Vr*Vr) + D/(5*Vr*Vr*Vr*Vr*Vr)
            + a13/(2*Tr*Tr*Tr*a15) * ( a14 + 1 - (a14+1+a15/(Vr*Vr))*exp(-a15/(Vr*Vr)) );

        phiCO2 = exp(lnphiCO2);

        return(phiCO2);
    }

    /**********************************************************************/
    /* Empirical model for pure water pressure (Duan 2003, Appendix B)    */
    /*                                                                    */
    /* this corresponds exactely to the result of the pwsat function      */
    /* in constrel2p2cni, but this function is faster                     */
    /**********************************************************************/

    double pure_water_pressure (double Temp)
    {
        double pcrit, Tcrit;
        double c1, c2, c3, c4, c5;
        double t, P;

        pcrit = 2.2085E7;     /* critical pressure [Pa] */
        Tcrit = 647.29;     /* critical temperature [K] */

        c1    = -38.640844;
        c2    = 5.8948420;
        c3    = 59.876516;
        c4    = 26.654627;
        c5    = 10.637097;

        t = (Temp - Tcrit)/Tcrit;

        P = pcrit*Temp/Tcrit * (1 + c1*pow(-t,1.9) + c2*t + c3*t*t + c4*t*t*t + c5*t*t*t*t);

        return P;
    }

};
}
#endif
