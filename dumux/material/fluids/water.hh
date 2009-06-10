// $Id$
#ifndef DUNE_FLUID_WATER_HH
#define DUNE_FLUID_WATER_HH

#include <dumux/material/property_baseclasses.hh>
#include <dumux/material/constrel/constrelwater.hh>

#include <dune/common/exceptions.hh>

namespace Dune
{
/*!
 * \ingroup properties
 *
 * \brief Fluid properties of water
 */
class Water : public Fluid
{
    ConstrelWater constRelWater;

public:
    Water(double constDensity = 0,
          double constViscosity = 0, double constEnthalpy = 0)
        : constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy)
    {}

    double viscosity (double T=283.15, double p=1e5, double X=1.0) const
    {
        if (constViscosity_)
            return constViscosity_;
        else
            return constRelWater.viscosity_water(T,p); //[kg/(ms)]
    }

    double density (double T=283.15, double p=1e5, double X=1.0) const
    {
        if (constDensity_)
            return constDensity_;
        else
            return 1000.0; // assumed to be incompressible[kg/m^3]
    }

    double enthalpy (double T=283.15, double p=1e5, double X = 1) const
    {
        if (constEnthalpy_)
            return constEnthalpy_;
        else {
            return constRelWater.enthalpy_water(T,p);
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

