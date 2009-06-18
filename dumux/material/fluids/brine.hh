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
#ifndef DUNE_FLUID_BRINE_HH
#define DUNE_FLUID_BRINE_HH

#include <dumux/material/property_baseclasses.hh>
#include <dumux/material/constrel/constrelbrine.hh>

#include <dune/common/exceptions.hh>

namespace Dune
{
/** \todo Please doc me! */

class Brine : public Fluid
{
    ConstrelBrine constRelBrine;

public:
    Brine(double constDensity = 0,
          double constViscosity = 0, double constEnthalpy = 0)
        :constDensity_(constDensity), constViscosity_(constViscosity), constEnthalpy_(constEnthalpy)
    {}

    double viscosity ( double T=283.15, double p=1e5, double X=1.) const
    {
        if (constViscosity_)
            return constViscosity_;
        else {
            double S;
            S = Salinity();
            return constRelBrine.viscosity_brine(T,S);
        }

        //           return 2.535e-4; // [kg/(ms)]
    }

    double Salinity() const
    {
        return 0.1;
    }

    double density ( double T=283.15, double p=1e5, double X=1.) const
    {
        if (constDensity_)
            return constDensity_;
        else {
            DUNE_THROW(NotImplemented, "Non-constant mass density of Brine");
        }
    }
    double enthalpy (double T=283.15, double p=1e5, double X = 1) const
    {
        if (constEnthalpy_)
            return constEnthalpy_;
        else {
            double S;
            S = Salinity();
            return constRelBrine.enthalpy_brine(T,p,S);
        }
    }

    double intEnergy(double T=283.15, double p=1e5, double X = 1) const
    {
        double intenergy;
        intenergy = enthalpy(T,p);
        return intenergy;
    }

private:

    double constDensity_;
    double constViscosity_;
    double constEnthalpy_;
};
}
#endif

