// $Id$

#ifndef DUNE_INTERSTITIALFLUID_TRAIL_HH
#define DUNE_INTERSTITIALFLUID_TRAIL_HH

#include <dumux/material/property_baseclasses.hh>

#include <dumux/material/constrel/constrelwater.hh>
#include <dune/common/exceptions.hh>

namespace Dune
{
/*!
 * \ingroup properties
 *
 * \brief Fluid properties of the Interstitial Fluid mixed with the
 *        TRAIL therapeutic agent.
 */
class InterstitialFluid_Trail : public Liquid_GL
{
public:
    InterstitialFluid_Trail()
    {}

    double viscosity (double T, double p, double X = 0.) const
    {
        return 0.00069152;
    }


    double density (double T, double p, double X = 0.) const
    {
        // TODO: correct?
        return 1025.0; // in [kg / m^3]
    }


    double enthalpy (double T, double p, double Xa = 0.) const
    {
        // TODO: correct enough?
        return ConstrelWater::enthalpy_water(T,p);
    }


    double intEnergy(double T, double p, double Xa = 0.) const
    {
        double u;
        double rho_mass = density(T,p);
        double h = enthalpy(T,p);

        u = h - (p / rho_mass);
        return u;
    }


    double diffCoeff(double T=283.15, double p=1e5) const
    {
        return 3.7378e-6;
    }


    double Xa_Max(double T, double p) const
    {
        return 0;
    }


    double p_vap(double T) const
    {
        return 0;
    }


    double henry(double T) const
    {
        return 0;
    }


    double T_vap(double p) const
    {
        return 0;
    }
};

}
#endif

