/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch, Andreas Lauser                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 * \brief Coarse configuration parameters for the lens problem.
 *
 * Just a few defines to allow quick experiments
 */
#ifndef STUPID_LENSCONFIG_HH
#define STUPID_LENSCONFIG_HH

// if defined, the boundary conditions and the material parameters are
// the same as for dumux' orignal lens problem.
//#define USE_ORIG_PROB

// if defined, most material parameters (capillary pressure, permeability,
// etc) are defined on the vertices of FE grid.
#define USE_VERTEX_PARAMETERS

// if defined, interface conditions of capillary pressure und and
// relative permeability are used at medium interfaces for cell based
// material parameters.
//#define USE_INTERFACE_CONDITION

// use parker-lenhard hysteresis for the simulation. if undefined only
// van-genuchten without hysteresis is used.
//#define USE_HYSTERESIS

#endif
