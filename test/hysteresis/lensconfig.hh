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
#ifndef DUMUX_LENSCONFIG_HH
#define DUMUX_LENSCONFIG_HH

// if defined, the boundary conditions and the material parameters are
// the same as for dumux' orignal lens problem.
#ifndef USE_ORIG_PROB
#define USE_ORIG_PROB 0
#endif

// if defined, most material parameters (capillary pressure, permeability,
// etc) are defined on the vertices of FE grid.
#ifndef USE_NODE_PARAMETERS
#define USE_NODE_PARAMETERS 1
#endif

// if defined, interface conditions of capillary pressure und and
// relative permeability are used at medium interfaces for element based
// material parameters.
#ifndef USE_INTERFACE_CONDITION
#define USE_INTERFACE_CONDITION 1
#endif

// use parker-lenhard hysteresis for the simulation. if undefined only
// van-genuchten without hysteresis is used.
#ifndef USE_HYSTERESIS
#define USE_HYSTERESIS 1
#endif

// use exponential splines to regularize the parker-lenhard hystersis model
#ifndef USE_SPLINES
#define USE_SPLINES 1
#endif

// use different parameters for the main imbibition and main drainage
// curves. if not set, the parker-lenhard model breaks down to plain
// van Genuchten with proper handling of the non-wetting phase
// residual saturation.
#ifndef USE_DIFFERENT_MAIN_CURVES
#define USE_DIFFERENT_MAIN_CURVES 0
#endif

// number of elements into X direction
#ifndef CELLRES_X
#define CELLRES_X 32
#endif

// number of elements into Y direction
#ifndef CELLRES_Y
#define CELLRES_Y 24
#endif

#endif
