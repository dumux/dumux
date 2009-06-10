// $Id$
#ifndef DUNE_FLUID_BRINE_HH
#define DUNE_FLUID_BRINE_HH

#include <dumux/material/property_baseclasses.hh>
#include <dumux/material/constrel/constrelbrine.hh>

#include <dune/common/exceptions.hh>

namespace Dune
{
/** \todo Please doc me! */

class Brine : public Fluid
{
    ConstrelBrine constRelBrine;

public:
    Brine(double constDensity = 0,
          double constViscosity = 0, double constEnthalpy = 0)
        :constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy)
    {}

    double viscosity ( double T=283.15, double p=1e5, double X=1.) const
    {
        if (constViscosity_)
            return constViscosity_;
        else {
            double S;
            S = Salinity();
            return constRelBrine.viscosity_brine(T,S);
        }

        //           return 2.535e-4; // [kg/(ms)]
    }

    double Salinity() const
    {
        return 0.1;
    }

    double density ( double T=283.15, double p=1e5, double X=1.) const
    {
        if (constDensity_)
            return constDensity_;
        else {
            double S, x_CO2_w;
            x_CO2_w = 0.0;
            S = Salinity();
            return constRelBrine.mass_density_brine_CO2(T,p,S,x_CO2_w);
        }
    }
    double enthalpy (double T=283.15, double p=1e5, double X = 1) const
    {
        if (constEnthalpy_)
            return constEnthalpy_;
        else {
            double S;
            S = Salinity();
            return constRelBrine.enthalpy_brine(T,p,S);
        }
    }

    double intEnergy(double T=283.15, double p=1e5, double X = 1) const
    {
        double intenergy;
        intenergy = enthalpy(T,p);
        return intenergy;
    }

private:

    double constDensity_;
    double constViscosity_;
    double constEnthalpy_;
};
}
#endif

