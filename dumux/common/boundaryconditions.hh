// $Id$
/*****************************************************************************
 *   Copyright (C) 2009-2010 by Bernd Flemisch                               *
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
#ifndef DUMUX_BOUNDARYCONDITIONS_HH
#define DUMUX_BOUNDARYCONDITIONS_HH

/**
* @file
* @brief  Definition of boundary condition types, extend if necessary
* @author Peter Bastian
*/
namespace Dumux
{
/** @addtogroup DISC_Operators
*
* @{
*/
/**
* @brief Define a class containing boundary condition flags
*
*/

//! base Class that defines boundary condition flags
struct BoundaryConditions
{
    /** \brief These values are ordered according to precedence */
    enum Flags {
        couplingOutflow = -2,
        couplingInflow = -1,
        outflow = 0,
        neumann = 1, //!< Neumann boundary
        process = 2, //!< Processor boundary
        dirichlet = 3 //!< Dirichlet boundary
    };
};

/** @} */
}
#endif
