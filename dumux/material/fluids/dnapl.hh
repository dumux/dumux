// $Id$
#ifndef DUNE_FLUID_DNAPL_HH
#define DUNE_FLUID_DNAPL_HH

#include <dumux/material/property_baseclasses.hh>
#include <dune/common/exceptions.hh>

namespace Dune
{
/** \ingroup properties
 * @brief Fluid properties of DNAPL
 * \todo which DNAPL?
 */
class DNAPL : public Fluid
{
public:
    DNAPL(double constDensity = 0,
          double constViscosity = 0, double constEnthalpy = 0)
        : constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy)
    {}

    double viscosity ( double T=283.15, double p=1e5, double X=1.) const
    {
        if (constViscosity_)
            return constViscosity_;
        else
            return 5.7e-4;//[kg/(ms)]
    }

    double density ( double T=283.15, double p=1e5, double X=1.) const
    {
        if (constDensity_)
            return constDensity_;
        else
            return 1460.0; // [kg/m^3]
    }
    double enthalpy (double T=283.15, double p=1e5, double X = 1) const
    {
        if (constEnthalpy_)
            return constEnthalpy_;
        else {
            //            return 1.0;
            DUNE_THROW(Dune::NotImplemented, "Non-constant enthalpy for DNAPL");
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

