// $Id$
#ifndef DUNE_FLUID_CO2_HH
#define DUNE_FLUID_CO2_HH

#include <dumux/material/property_baseclasses.hh>
#include <dumux/material/constrel/constrelco2.hh>

#include <dune/common/exceptions.hh>

namespace Dune
{
/** \ingroup properties
 * @brief Fluid properties of CO2
 */
class CO2 : public Fluid
{
    ConstrelCO2 constRelCO2;

public:
    CO2(double constDensity = 0,
        double constViscosity = 0, double constEnthalpy = 0)
        : constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy)
    {}

    double density ( double T=283.15, double p=1e5, double X=1.) const
    {
        if (constDensity_)
            return constDensity_;
        else
            return constRelCO2.density(T,p);
    }
    double viscosity ( double T=283.15, double p=1e5, double X=1.) const
    {
        if (constViscosity_)
            return constViscosity_;
        else
            return constRelCO2.viscosity(T,p);
    }

    double enthalpy ( double T=283.15, double p=1e5, double X = 1) const
    {
        if (constEnthalpy_)
            return constEnthalpy_;
        else {
            return constRelCO2.enthalpy(T,p);
        }
    }

    double intEnergy( double T=432, double p=3.086e7, double X = 1) const
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

