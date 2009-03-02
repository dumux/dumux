// $Id$

#ifndef DUNE_CONSTRELCO2_HH
#define DUNE_CONSTRELCO2_HH

#include "co2tables.hh"

#include <dumux/exceptions.hh>

#include <math.h>

namespace Dune
{

/** \todo Please doc me! */

class ConstrelCO2
{
public:
    // from MUFTE
    // calcualted from Span and Wagner (1996)
    static double density(double T, double p)
    {
        return Dune::Co2Tables::tabulatedDensity.at(p, T);
    }


    // from MUFTE:
    /*******************************************************************/
    /*                                                                 */
    /* Viscosity of CO2: - Vesovic et al., 1990                        */
    /*                   - Fenghour et al., 1998                       */
    /*                                                                 */
    /*******************************************************************/
    static double viscosity(double temp, double pg)
    {
        static const double a0 = 0.235156;
        static const double a1 = -0.491266;
        static const double a2 = 5.211155E-2;
        static const double a3 = 5.347906E-2;
        static const double a4 = -1.537102E-2;

        static const double d11 = 0.4071119E-2;
        static const double d21 = 0.7198037E-4;
        static const double d64 = 0.2411697E-16;
        static const double d81 = 0.2971072E-22;
        static const double d82 = -0.1627888E-22;

        static const double ESP = 251.196;

        double mu0, SigmaStar, TStar;
        double dmu, rho;
        double visco_CO2;

        if(temp < 275.) // regularisation
        {
            // temp = 275;
            DUNE_THROW(Dune::NumericalProblem,
                       "ConstrelCO2: Temperature " << temp << " out of range at " << __FILE__ << ":" << __LINE__);
        }


        TStar = temp/ESP;

        /* mu0: viscosity in zero-density limit */
        SigmaStar = exp(a0 + a1*log(TStar)
                        + a2*log(TStar)*log(TStar)
                        + a3*log(TStar)*log(TStar)*log(TStar)
                        + a4*log(TStar)*log(TStar)*log(TStar)*log(TStar) );

        mu0 = 1.00697*sqrt(temp) / SigmaStar;

        /* dmu : excess viscosity at elevated density */

        rho = density(temp, pg); /* CO2 mass density [kg/m^3] */

        dmu = d11*rho + d21*rho*rho + d64*pow(rho,6)/(TStar*TStar*TStar)
            + d81*pow(rho,8) + d82*pow(rho,8)/TStar;

        /* dmucrit : viscosity increase near the critical point */

        // False (Lybke 2July2007)
        //e1 = 5.5930E-3;
        //e2 = 6.1757E-5;
        //e4 = 2.6430E-11;
        //dmucrit = e1*rho + e2*rho*rho + e4*rho*rho*rho;
        //visco_CO2 = (mu0 + dmu + dmucrit)/1.0E6;   /* conversion to [Pa s] */

        visco_CO2 = (mu0 + dmu)/1.0E6;   /* conversion to [Pa s] */

        return visco_CO2;
    }

    static double enthalpy(double temp, double pg)
    {
        return Dune::Co2Tables::tabulatedEnthalpy.at(pg, temp);
    }

    // heat conductivity
    static double lambda(double temp, double pg)
    {
        DUNE_THROW(NotImplemented, "Heat conductivity of CO2");
    }

    static double satpressure(double temp)
    {
        /*values from Span und Flacke 2004: Statusbericht DKV Nr 20*/
        static const double ts[] = {
            270, 272, 276, 280, 284,
            288, 292, 296, 300, 304
        };
        static const double ps[] = {
            3.2033E6, 3.3802E6, 3.7555E6, 4.1607E6, 4.5978E6,
            5.0688E6, 5.5761E6, 6.1227E6, 6.7131E6, 7.3555E6
        };

        /*                if(temp < 270 || temp > 304)
                          {
                          DUNE_THROW(Dune::NumericalProblem,
                          "ConstrelCO2: Temperature " << temp << " out of range at " << __FILE__ << ":" << __LINE__);
                          }
        */

        if (temp > 304.0)
            temp = 304.0;

        for (int i=0; i < 9; i++)
        {
            if (ts[i] <= temp && temp <= ts[i+1])
            {
                return ps[i] +
                    (temp - ts[i])/(ts[i+1] - ts[i])
                    * (ps[i+1] - ps[i]);
            }
        }

        DUNE_THROW(Dune::NumericalProblem, "satpressure_CO2: Temperature " << temp << " out of range");
    }
};
};

#endif
