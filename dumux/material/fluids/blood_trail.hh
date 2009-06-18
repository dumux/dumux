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
    Blood_Trail()
    {}

    double viscosity (double T, double p, double X = 0.) const
    {
        return 0.0069152;
    }


    double density (double T, double p, double X = 0.) const
    {
        // TODO: correct?
        return 1; // in [kg/m^3]
    }


    double enthalpy (double T, double p, double Xa = 0.) const
    {
    	DUNE_THROW(NotImplemented,
    	                   "Enthalpy of TRAIL in human blood");
    }


    double intEnergy(double T, double p, double Xa = 0.) const
    {
    	DUNE_THROW(NotImplemented,
    	                   "Internal energy of TRAIL in human blood");
    }


    double diffCoeff(double T=283.15, double p=1e5) const
    {
        return 3.738e-6;
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

