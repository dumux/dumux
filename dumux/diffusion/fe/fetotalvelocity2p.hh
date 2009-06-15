#ifndef DUNE_FETOTALVELOCITY2P_HH
#define DUNE_FETOTALVELOCITY2P_HH
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
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
namespace Dune
{

template<class GridView, class Scalar, class VC, class Problem, class LocalStiffnessType, class Communication>
void FEPressure2PBase<GridView, Scalar, VC, Problem, LocalStiffnessType, Communication>::calculateVelocity(const Scalar t=0) const
{
    DUNE_THROW(Dune::NotImplemented, "velocities only implemented for finite volume and mimetic finite differences discretisations");
    return;
}

}
#endif
