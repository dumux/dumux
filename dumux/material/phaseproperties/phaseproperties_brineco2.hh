/// $Id$

#ifndef PHASEPROPERTIES_BRINECO2_HH_
#define PHASEPROPERTIES_BRINECO2_HH_

#include <dumux/material/property_baseclasses.hh>
#include <dumux/material/constrel/constrelbrine.hh>
#include <dumux/material/constrel/constrelco2.hh>
#include <dumux/material/constrel/constrelair.hh>
#include <math.h>

/**
 * \ingroup properties
 * \author Jochen Fritz
 */

namespace Dune
{

/**\ingroup properties
 * @brief property class for gaseous phase of air/water mixture
 */
class Gas_BrineCO2 : public Gas_GL
{
public:
	virtual double density(double T, double p, double Xw=0.) const // [kg / m^3]
	{
		return constRelCO2.density(T,p); // see constrelco2.hh
	}

	virtual double viscosity(double T, double p, double Xw=0.) const
	{
        return constRelCO2.viscosity(T, p); // see constrelco2.hh
	}

/*	virtual double viscosityCO2(double T, double p, double rho, double Xw=0.) const // [kg / (m*s)]
	{
		return constRelCO2.viscosity(T,p,rho); // see constrelco2.hh
	}
*/

	virtual double intEnergy(double T, double p, double Xw=0.) const
	{
    	double u;
    	double rho_mass = density(T,p);
    	double h = enthalpy(T,p);

    	u = h - (p / rho_mass);
    	return u;
	}

	virtual double enthalpy(double T, double p, double Xw=0.) const
	{
		return constRelCO2.enthalpy(T,p);
	}

	virtual double diffCoeff(double T=283.15, double p=1e5) const
	{
		return 2.6e-5; // [m^2/s] diffusion coefficient of water in gas phase
	}

	virtual double Xw_Max(double T, double p) const
	{
		return 0.001;
	}

	Gas_BrineCO2() : Gas_GL()
		{}
private:
	ConstrelCO2 constRelCO2;
};

/**\ingroup properties
 * @brief property class for liquid phase of air/water mixture
 */
class Liq_BrineCO2 : public Liquid_GL
{
public:
	virtual double density(double T, double p, double Xa=0.) const // [kg / m^3]
	{
		return constRelBrine.mass_density_brine_CO2(T,p,Salinity,Xa);
	}

	virtual double viscosity(double T, double p, double Xa=0.) const
	{
		return constRelBrine.viscosity_brine(T,Salinity);
	}

	virtual double intEnergy(double T, double p, double Xa=0.) const
	{
		return constRelBrine.enthalpy_brine(T,p,Salinity);
	}

	virtual double enthalpy(double T, double p, double Xa=0.) const
	{
		return constRelBrine.enthalpy_brine(T,p,Salinity);
	}

	virtual double diffCoeff(double T=283.15, double p=1e5) const
	{
		return 2.E-9; // [m^2/s] diffusion coefficient of CO2 in liquid phase

	}

	virtual double henry(double T) const
	{

		return 0.; // [1/Pa]
	}

	virtual double Xa_Max(double T, double p) const
	{
		return 0.;
	}

	virtual double p_vap(double T) const
	{
		return 0.; //[Pa]
	}

	virtual double T_vap(double p) const
	{
		return 0.;
	}

  Liq_BrineCO2() : Liquid_GL()
	{
  	Salinity = 0.1;
	}

private:
	ConstrelBrine constRelBrine;
	ConstrelAir constRelAir;
	double Salinity;


};

} // namespace

#endif /*PHASEPROPERTIES_HH_*/
