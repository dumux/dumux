/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
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
 *
 * \brief Performs integration using implicit euler.
 */
#ifndef DUMUX_IMPLICITEULERSTEP_HH
#define DUMUX_IMPLICITEULERSTEP_HH

#include <algorithm>

namespace Dune {
    /*!
     * \brief Performs integration using implicit euler.
     */
    template<class Problem>
    class ImplicitEulerStep
    {
        typedef typename Problem::DomainTraits::Scalar Scalar;

    public:
        //! excute an implicit euler integration.
        void execute(Problem &problem,
                     Scalar t,
                     Scalar &dt,
                     Scalar &nextDt,
                     Scalar maxDt,
                     Scalar tEnd,
                     Scalar cFLFactor)
            {
                Scalar eps = 1e-8;
                dt = std::min( dt, maxDt );
                if (tEnd - t <= (1+eps)*dt)
                    dt = (1+eps)*(tEnd - t);

                problem.updateModel(dt, nextDt);

                nextDt = std::min(nextDt, tEnd - t);
                nextDt = std::min(nextDt, maxDt);
            }
    };
}
#endif
