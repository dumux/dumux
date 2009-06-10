// $Id$

#ifndef DUNE_BLOOD_TRAIL_HH
#define DUNE_BLOOD_TRAIL_HH

#include <dumux/material/property_baseclasses.hh>
#include <dune/common/exceptions.hh>

namespace Dune
{
/*!
 * \ingroup properties
 *
 * \brief Fluid properties of human blood mixed with the TRAIL
 *        therapeutic agent.
 */
class Blood_Trail : public Liquid_GL
{
public:
    Blood()
    {}

    double viscosity (double T, double p, double X = 0.) const
    {
        return 0.0069152;
    }


    double density (double T, double p, double X = 0.) const
    {
        // TODO: correct?
        return 1060.0; // in [kg/m^3]
    }


    double enthalpy (double T, double p, double Xa = 0.) const
    {
        // TODO: correct enough?
        return constRelWater.enthalpy_water(T,p);
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
        return 2.7378e-6;
    }


    double Xa_Max(double T, double p) const
    {
        DUNE_THROW(NotImplemented, 
                   "Maximum mass concentration of TRAIL in human blood");
    }


    double p_vap(double T) const
    {
        DUNE_THROW(NotImplemented, 
                   "Vapour pressure of TRAIL in human blood");
    }


    double henry(double T) const
    {
        DUNE_THROW(NotImplemented, 
                   "Henry constant of TRAIL in human blood");
    }


    double T_vap(double p) const
    {
        DUNE_THROW(NotImplemented, 
                   "Vapourization temperature of TRAIL in human blood");
    }
};

}
#endif

