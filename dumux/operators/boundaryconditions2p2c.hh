// $Id$
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

#ifndef DUNE_BOUNDARYCONDITIONS2P2C_HH
#define DUNE_BOUNDARYCONDITIONS2P2C_HH


namespace Dune
{
/** \ingroup decoupled2p2c
 * \brief Defines type of boundary conditions for 2p2c processes
 *
 *  This is to distinguish BC types for 2p2c processes similar to
 * the class Dune::BoundaryConditions which distinguishes between
 * dirichlet, process and neumann.
 * BoundaryConditions for compositional transport can either be
 * defined as saturations or as total concentrations (decoupled method).
 * Each leads via pressure and constitutive relationships to the other
 * and to the phase concentrations.
 */
struct BoundaryConditions2p2c
{
    enum Flags
        {
            saturation=1,
            concentration=2,
        };
};

}

#endif
