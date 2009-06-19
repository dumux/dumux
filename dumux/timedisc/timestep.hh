// $Id$
/*****************************************************************************
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
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

#ifndef DUNE_TIMESTEP_HH
#define DUNE_TIMESTEP_HH

namespace Dune {

/** \todo Please doc me! */

template<class G, class Model>
class TimeStep
{
public:
    virtual void execute(Model& model, double t, double& dt,
                         double maxDt, double tEnd, double cFLFactor) = 0;

    //! always define virtual destructor in abstract base class
    virtual ~TimeStep () {}
};
}
#endif
