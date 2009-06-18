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

#ifndef DUNE_IMPLICITEULERSTEP_HH
#define DUNE_IMPLICITEULERSTEP_HH

namespace Dune {

/** \todo Please doc me! */

template<class G, class Model, bool useOldNewton = true>
class ImplicitEulerStep : public TimeStep<G, Model>
{
public:
    void execute(Model& model, double t, double& dt,
                 double maxDt, double tEnd, double cFLFactor)
    {
        double eps = 1e-8;
        dt = std::min( dt, maxDt );

        if (tEnd - t <= (1+eps)*dt)
            dt = (1+eps)*(tEnd - t);

        model.update(dt);
        dt = std::min(dt, tEnd - t);

        return;
    }

};

template<class G, class Model>
class ImplicitEulerStep<G, Model, false> : public TimeStep<G, Model>
{
public:
    void execute(Model& model, double t, double& dt,
                 double maxDt, double tEnd, double cFLFactor)
    {
        double eps = 1e-8;
        dt = std::min( dt, maxDt );

        if (tEnd - t <= (1+eps)*dt)
            dt = (1+eps)*(tEnd - t);

        double nextDt;
        model.updateModel(dt, nextDt);
        dt = std::min(nextDt, tEnd - t);

        return;
    }

};
}
#endif



