// $Id$
#ifndef DUNE_PUREFLUIDS_HH
#define DUNE_PUREFLUIDS_HH

#include <dumux/material/property_baseclasses.hh>
#include <dumux/material/constrel/constrelco2.hh>
#include <dumux/material/constrel/constrelwater.hh>
#include <dumux/material/constrel/constrelbrine.hh>
#include <dumux/material/constrel/constrelair.hh>

#include <dune/common/exceptions.hh>

namespace Dune
{
/** \ingroup properties
 * @brief Fluid properties of methane
 *
 * \todo this class uses the constant relations of water too much! it
 *       is definitely wrong!
 */
class Methane : public Fluid
{
    ConstrelWater constRelWater;

public:
    Methane(double constDensity = 0,
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
        {
            double R = 0.54978284;
            double temp = 284.43;
            double rhog = p/(R*temp*1000.0);
            if(rhog < 1e-15)
                rhog = 1e-15;

            return(rhog);
        }
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

