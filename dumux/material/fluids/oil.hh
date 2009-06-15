// $Id$
#ifndef DUNE_FLUID_OIL_HH
#define DUNE_FLUID_OIL_HH

#include <dumux/material/property_baseclasses.hh>

#include <dune/common/exceptions.hh>

namespace Dune
{
/** \todo Please doc me! */

class Oil : public Fluid
{
public:
    Oil(double constDensity = 0,
        double constViscosity = 0, double constEnthalpy = 0)
        : constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy)
    {}

    double viscosity ( double T=283.15, double p=1e5, double X=1.0) const
    {
        if (constViscosity_)
            return constViscosity_;
        else
            return 800e-3;//[kg/(ms)]
    }

    double density ( double T=283.15, double p=1e5, double X=1.0) const
    {
        if (constDensity_)
            return constDensity_;
        else
            return 890.0; // [kg/m^3]
    }
    double enthalpy (double T=283.15, double p=1e5, double X = 1) const
    {
        if (constEnthalpy_)
            return constEnthalpy_;
        else {
            //            return constRelOil.enthalpy(T,p,X);
            // TODO
            DUNE_THROW(Dune::NotImplemented, "Non-constant enthalpy of oil");
        }
    }

    double intEnergy( double T=283.15, double p=1e5, double X = 1) const
    {
        double u;
        double rho_mass = density(T,p);
        double h = enthalpy(T,p);

        u = h - (p / rho_mass);
        return u;
    }

private:
    double constDensity_;
    double constViscosity_;
    double constEnthalpy_;
};
}

#endif

