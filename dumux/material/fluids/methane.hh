// $Id:$
/*****************************************************************************
 *   Copyright (C) <YEARS> by <ADD_AUTHOR_HERE>                              *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUNE_PUREFLUIDS_HH
#define DUNE_PUREFLUIDS_HH

#include <dumux/material/property_baseclasses.hh>
#include <dumux/material/constrel/constrelair.hh>

#include <dune/common/exceptions.hh>

namespace Dune
{
/** \ingroup properties
 * @brief Fluid properties of methane
 *
 * \todo this class uses the constant relations of air quite a lot! 
 *       is it wrong?
 */
class Methane : public Fluid
{
    ConstrelAir constRelAir;

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
            return constRelAir.viscosity_air(T,p); //[kg/(ms)]
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
            // \todo 
            DUNE_THROW(NotImplemented, "Enthalpy of air");
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

