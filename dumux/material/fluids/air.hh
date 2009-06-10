// $Id$
#ifndef DUNE_FLUID_AIR_HH
#define DUNE_FLUID_AIR_HH

#include <dumux/material/property_baseclasses.hh>
#include <dumux/material/constrel/constrelair.hh>

#include <dune/common/exceptions.hh>

namespace Dune
{
/** \ingroup properties
 * @brief Fluid properties of Air
 */
class Air : public Fluid
{
    ConstrelAir constRelAir;

public:
    Air(double constDensity = 0,
        double constViscosity = 0, double constEnthalpy = 0)
        :constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy)
    {}

    double viscosity ( double T=283.15, double p=1e5, double X=0.0) const
    {
        if (constViscosity_)
            return constViscosity_;
        else
            return constRelAir.viscosity_air(T); //[kg/(ms)]
    }

    double density ( double T=283.15, double p=1e5, double X=0.0) const
    {
        if (constDensity_)
            return constDensity_;
        else {
            const double molarMassAir = 0.02896; // [kg/mole]

            return constRelAir.rho_idGG_mass(T,p,molarMassAir); // [kg/m^3]
        }
    }
    double enthalpy (double T=283.15, double p=1e5, double X = 0.0) const
    {
        if (constEnthalpy_)
            return constEnthalpy_;
        else {
            return constRelAir.sp_enth2p2cni_g(T,p,0.0);
        }
    }

    double intEnergy( double T=283.15, double p=1e5, double X = 0.0) const
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

