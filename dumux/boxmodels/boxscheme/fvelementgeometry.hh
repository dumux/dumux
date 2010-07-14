// $Id: fvelementgeometry.hh 3774 2010-06-24 06:22:44Z bernd $
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
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

#ifndef DUMUX_FVELEMENTGEOMETRY_HH
#define DUMUX_FVELEMENTGEOMETRY_HH

#if HAVE_DUNE_PDELAB
#include "fvelementgeometry-pdelab.hh"
#else
#include "fvelementgeometry-disc.hh"
#endif

#endif

